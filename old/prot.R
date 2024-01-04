
## dir <- "~/pap/21/giovanni-exch/prot/R/" # just for P
## setwd(dir)
## source("~/s/first.s")

rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.15
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(135)

library(dplyr)
library(tidyr)
library(mvtnorm)
library(invgamma)
library(splines)
## library(clara)
## create data as
## y:  (R x C) with R=32 rows and C=3990 columns
## X:  (R x p): design matrix
## R=32   patients
## C=3990 proteins
##

readDta = function(){
    ## cleaning data and construct design matrix
    ## rownames(dta)=c()
    load("../Data-and-Results/data_protein.RData")
    dta=cbind(c(rep(0,21),rep(1,21)),PL)
    colnames(dta)[1]="z"

    dta2=as.data.frame(dta) %>%
        group_by(z,age) %>%
        summarise(across(everything(), mean, na.rm = TRUE)) %>%
        as.matrix()
    ages <- matrix(dta2[,2],ncol=2)
    ## remove proteins with < 10 not NA's
    nna  <- apply(dta2,2,function(x){sum(is.na(x))})
    R = nrow(dta2)
    jdx  <- which (30-nna<10)
    dta2=dta2[,-jdx]
    nna  <- apply(dta2,1,function(x){sum(is.na(x))})
    ## could drop obs with too many NA's -- but don't do it now..
    
    age=dta2[,2]
    age_std=(dta2[,2]-mean(dta2[,2]))/sd(dta2[,2])
    x=cbind(dta2[,1],age_std,dta2[,1]*age_std)
    ## don't think we ever use the last item, the interaction
    y= dta2[,3:dim(dta2)[2]]
    ## remove some more proteins with too many NA's;
    ## will use hclust(.) to initialize clustering; therefore remove few more proteins that
    ## give NA's in evaluating dist(.) (needed for hclust):
    idx  = unique(which(is.na(as.matrix(dist(t(y)))),arr.ind=T)[,1]) #2nd col is the same..
    y = y[ ,-idx]
    ## tried the same for patients --
    ##        but turns out hclust is ok with all patients as is
    
    ## total number of proteins
    C=dim(y)[2]

    ## construct splines using bs package
    bsx2=bs(x[,2],degree = 3,
            knots = quantile(x[,2],c(0.33,0.67)),intercept = T)
    X=cbind(bsx2,x[,1]*bsx2)

    return(list(X=X,y=y,C=C,R=nrow(X),p=ncol(X), ages=ages))
}

## global variables for the dta
out = readDta() ## make global vars for the data
X = out$X
y = out$y
my = mean(y,na.rm=T)
sy =  mean(apply(y,2,var,na.rm=T))
ages = out$ages ## age grid (for plotting)
## boxplot(apply(y,2,mean,na.rm=T))
## boxplot(sqrt(apply(y,2,var,na.rm=T)))

C=out$C
R=out$R
p=out$p
## R= # patients (=# time points x2); C=  #proteins


## prior pars
prior = list(
    ## column effects (proteins)
    meta = rep(0,p),
    Seta = diag(p),  ## might want to reduce var of trt offsets (2nd part)
    aeta = 1,        ## eta[j] ~ G_eta; G_eta ~ DP(aeta, N(meta,Seta))
    beta = 0.05,
    ## row effects (patients)
    mxi = my,       ## xi[i] ~ G_xi; G_xi ~ DP(axi, N(mxi,Sxi))
    Sxi = 1,
    axi = 1,
    bxi = 0.05,
    ## hyperpars theta=sig2
    asig = 1,
    bsig = 5) ## w=1/sigs ~ Ga(asig/2, bsig/2); E(w) = m=a/b, V(w)=m/(b/2)


    
    

main = function(niter=100){
    ## main driver
    sigs = 1/ (prior$asig/prior$bsig) ## prior mean for residual var
    mdpEta = mdpInit(2,C,q=p, mbeta=prior$meta, Sbeta=prior$Seta,
                     a=prior$aeta, b=prior$beta, sigs=sigs)
    mdpXi =  mdpInit(1,R,q=1, mbeta=prior$mxi, Sbeta=prior$Sxi,
                     a=prior$axi, b=prior$bxi, sigs=sigs)
    ## initialize partitions using hclust;
    ## in mdpInit they were initialized with K=1 cluster
    mdpOffset(mdpEta)
    mdpXi  = mdpInitClust(mdpXi, R) # all singleton clusters
    mdpOffset(mdpXi)
    Kmax = if(mdpEta$b>0) round( mdpEta$a/mdpEta$b ) else 20
    mdpEta = mdpInitClust(mdpEta, Kmax)
    
    mcmcInit()
    for(it in 1:niter){
        mdpOffset(mdpXi)  # updates global yt = y - pat effects
        mdpEta = mdpUpdate(mdpEta,niter=500)
        mdpOffset(mdpEta)  # updates global yt = y - prot effects
        mdpXi = mdpUpdate(mdpXi)
        mdpOffset(mdpXi,incr=T)  # yt = y - (pat effects + prot effects)
        mdpXi$sigs = updateSigs(mdpXi)
        mdpEta$sigs = mdpXi$sigs    # need residual var in both
        if (it %% 10==0)
            mcmcUpd(mdpEta,mdpXi,it)
    }# it
}

mcmcInit = function(){
    ## global variables for debugging and summaries
    Ey <<- matrix(0,nrow=R,ncol=C)     # running update for fitted values
    Ey2 <<- Ey      # .. for 2nd moment (for sd)
    Eyp <<- Ey
    Eyp2 <<- Ey
    nEy <<- 0        # number of updates (Ey and Ey2 are running sums -- divide by nEy for plotting)
    SSM <<- 0        # mean residual SS
    chain <<- NULL # collect summaries at each iteration
    ## it, SSM, sig2, K-pat, K-prot, nk-pat*, nk-prot* (nk*: 5 largest cluster sizes, padded with 0's if needed)
    filesAppend <<- F # init the files at first write
    return(0)
}
    
mcmcUpd = function(mdpEta, mdpXi,it){
    ## updates global summaries and writes to file (if desired..)
    ## prints out summary at iteration

    ## print summaries
    cat(it, "\t K-Prot=",mdpEta$K," (",table(mdpEta$s),"), \n avg(betas)=",
        format(apply(mdpEta$betas,1,mean),digits=2),"\n")
    cat("\t K-Pat= ",mdpXi$K," (",table(mdpXi$s),") avg(betas)=",
        format(mean(mdpXi$betas),digits=2),"\t")
    cat("sig=", format(sqrt(mdpXi$sigs)),"   SSM=", format(SSM), "\n")

    ## summaries
    nkProt = sort( table(mdpEta$s), dec=T)[1:5]
    nkPat  = sort( table(mdpXi$s),  dec=T)[1:5]
    nkProt = ifelse(is.na(nkProt),0,nkProt) # replace NA's by 0
    nkPat = ifelse(is.na(nkPat),0,nkPat)

    line = c(it, SSM, mdpXi$sigs, mdpEta$K, mdpXi$K, nkProt, nkPat)
    names(line) = c("it", "SSM", "sig2", "K-pat", "K-prot", paste("nk",1:5,sep=""), paste("nk",1:5,sep=""))
    chain <<- rbind(chain,line)
    ## it, SSM, sig2, K-pat, K-prot, nk-pat*, nk-prot* (nk*: 5 largest cluster sizes, padded with 0's if needed)

    ## updated fitted summaries
    yhat = mdpFitted(mdpEta, mdpXi)
    yphat = mdpFitted(mdpEta, fitProt=T)
    Ey <<- Ey+yhat
    Ey2 <<- Ey2+yhat*yhat
    Eyp <<- Eyp+yphat
    Eyp2 <<- Eyp2+yphat*yphat
    nEy <<- nEy+1
    browser()
    if (it %% 50 == 0){
        options(digits=2)
        write.table(chain,"chain.txt",sep=",",append=filesAppend)
        write.table(mdpXi$s,"sPat.txt",sep=",",append=filesAppend)
        write.table(mdpEta$s[1:250],"sProt.txt",sep=",",append=filesAppend)
        ## save only the first 250 proteins -- seems enough :-)
        filesAppend <<- TRUE # from now on append those three files..
                             ## the files below are always overwritten
        write.table(format(Ey/nEy),"Ey.txt",quote=F,col.names=F,row.names=F,sep=",")
        write.table(format(Ey2/nEy),"Ey2.txt",quote=F,col.names=F,row.names=F,sep=",")
        write.table(format(Eyp/nEy),"Eyp.txt",quote=F,col.names=F,row.names=F,sep=",")
        write.table(format(Eyp2/nEy),"Eyp2.txt",quote=F,col.names=F,row.names=F,sep=",")
    }
    return(0)
}

mdpFitted = function(mdpEta, mdpXi=NULL, fitProt=F){
    ## fitProt: indicator for fitting *without* patient effects
    f = matrix(0, nrow=R, ncol=C)
    fitPat = 0 # initialize, vector of patient-specific means
    if (!fitProt)
        fitPat = mdpXi$betas[,mdpXi$s] # patient-specific mean
    for(k in 1:mdpEta$K){
        idx = which(mdpEta$s==k)
        f[,idx] = fitPat + X%*%mdpEta$betas[,k] # add prot-spec spline
    }
    return(f)
}
    


mdpOffset = function(mdp,incr=F){
    ## creates global variable
    ## yt = y-xxx, where xxx are either protein effects (mdp$d=2),
    ##   or patient effects (d=1), or both (mdp2 != NULL)
    ##   (the latter is used to update residual var
    ## if correcting for both, then mdp2 must be mdpXi!
    if (mdp$d==1){ # correcting for patient effects
        offset = mdp$betas[,mdp$s]
    } else { # d=2, correcting for protein effects
        offset = 0*y         # initialize matrix of right size
        for(k in 1:mdp$K){
            idx = which(mdp$s==k)
            offset[,idx] = X%*%mdp$betas[,k]
        }
    } # correct for protein effects (d=2)
    if (incr)            ## additional correction (use for residual var)
        yt <<- yt-offset ## additional correction
    else 
        yt <<- y-offset  ## under d=1: repeating offset for each column (=prot)
    return(0)
}


#######################################################
## DP mix of normal linear regressions

## use list with
## global vars for data:
## y:  (RxC) R=# patients, C=# proteins
## offset:  (R x C)
## yt: (RxC) copy of y, for (y-offset) 
## X:  (Rxp)
##
## parameters:
## d:       dimension: d=1 for patients (="rows"); d=2 for proteins (=cols)
## M:       # experimental units (under d=1: R; under d=2: C)
## q:       dimension of cluster-specific pars (1 under d=1; p under d=2)
## a,b:     hyperpar's for PY (b=0 for DP)
## mbeta,Sbeta:   (qx1) and (qxq) mv or univ normal prior moments for etas[k]
## sigs:    residual variance
##
## local variables
## s:       (Mx1) cluster membership indicators
## K:       # clusters
## betas:   (q x K) cluster specific pars (columns= cluster)

## global variables (for debugging and summaries):
## MSS:     mean residual
## Ey
## Ey2
## nupd

mdpInit = function(d,M,q,
                   mbeta=rep(0,q),
                   Sbeta=diag(q),
                   a=1,b=0,
                   sigs=1){
    ## initializes pars for a DP mix of normal lin regression
    ## use initially K=1 cluster
    Lbeta = solve(t(chol(Sbeta))) # Leta*S*Leta' = I :-)
    Sibeta = solve(Sbeta)
    Simbeta = Sibeta%*%mbeta
    mdp = list(d=d,M=M,q=q,
               s=rep(1,M), K=1,a=a,b=b,
               betas=as.matrix(mbeta,ncol=1),
               mbeta=mbeta, Sbeta=Sbeta, Sibeta=Sibeta,
               Lbeta=Lbeta, Simbeta=Simbeta,
               sigs=sigs)
    return(mdp)
}

mdpSwap = function(mdp,k1,k2){
    ## swaps clusters k1 vs k2
    if (k1 != k2){
        idx1 = which(mdp$s==k1)
        idx2 = which(mdp$s==k2)
        mdp$s[idx1] = k2
        mdp$s[idx2] = k1
        b1 = mdp$betas[,k1]
        mdp$betas[,k1] = mdp$betas[,k2]
        mdp$betas[,k2] = b1
    }
    return(mdp)
}
    
mdpInitClust = function(mdp,K0){
    ## initialize clustering using clusters from hclust
    if (mdp$d==1){ # clustering patients, ignore K0 :-)
        mdp$s=1:R
        mdp$K=R
    } else{ 
        D = dist(t(y))
        hc = hclust(D)
        mdp$s = cutree(hc,k=K0)
        mdp$K = length(unique(mdp$s))
    }
    mdp = mdpUpdB(mdp,spl=F)
    return(mdp)
}

mdpUpdate = function(mdp, niter=mdp$M){
    mdp = mdpUpdS(mdp,niter)
    mdp = mdpUpdB(mdp)
    return(mdp)
}

mdpUpdS = function(mdp, niter=mdp$M, ns=20){
    ## MCMC transition probs to update s (if updateS=T);
    ## niter iterations
    ## 1. evaluate fk=p(y*[k] | s) (marginalizing eta), k=1..K
    ## 2. iteratively update one random s, updating fk as needed

    ## 1. evaluate fk=p(y*[k] | s) (marginalizing eta), k=1..K
    ## to evaluate likel of (currently considered) unit i in cluster k we
    ## 1a. for n[k]<ns we use marg likel p(y[i] | y*[k]-)=p(y*[k])/p(y*[k]-)
    ##     here marginalization is w.r.t. betas[k],
    ##     and y*[k]- = cluster k w/o unit i ("llwo")
    ##         y*[k]  = .. with unit i ("llwith")
    ## 1b. otherwise conditional             p(y[i] | betas[k]) (much faster)

    ## 1. evaluate marginal (only for nk[k]<ns)
    llwith = rep(0,mdp$K)   # log cluster-spec likel *with* (currently considered) unit i
    llwo   = rep(0,mdp$K)   # .. *without*
    nk = rep(0,mdp$K)       # cluster sizes
    for(k in 1:mdp$K){ 
        nk[k] = sum(mdp$s==k)
        if (nk[k]<ns)
            llwo[k] = mdpMarg(k,mdp) # initialize llwo
        ## note for k=s[i] llwo[k] is actually ll"with"! Will correct this below (*)
    }
    llwith[nk >= ns] = NA     ## do *not* evaluate llwith & llwo for n[k]>ns
    llwo[nk>=ns] = NA
    ## 2. iteratively update one random s, updating fk as needed
    iit = sample(1:mdp$M,size=niter,replace=T)
       ## randomly chosen unit i to update s[i]
    for(it in 1:niter){
        i = iit[it]
        ## 2.1 take i out of it's current cluster
        sold = mdp$s[i]        # save current s[i]
        mdp$s[i]= -1           # take it out of current cluster
        nki = nk[sold]         # current cluster size for i
        nk[sold]  = nki-1
        if (nk[sold]==0){      # singleton cluster -- take it out and decrement K
            mdp = mdpSwap(mdp,sold,mdp$K)
            llwith[mdp$K] = llwo[sold] # currently not used - but save it
            llwo[sold] = llwo[mdp$K]
            nk[sold] = nk[mdp$K]
            nk = nk[-mdp$K]            # drop last element
            mdp$K = mdp$K-1
        } else {
            if (nk[sold]<ns){      # llwo was not evaluated; get it..
                llwith[sold] = llwo[sold] # currently not used - but save it
                llwo[sold] = mdpMarg(sold,mdp) # evaluate (new) llwo (w/o i)
            }
        }
        lps = rep(0,mdp$K)     # initialize (log) post prob p(s[i]=k | ...)
        prnew = mdp$a-mdp$b*mdp$K
        lps       = log(nk+mdp$b)
        if (prnew > 0){
            lps[mdp$K+1] = log(prnew) ## last element for new cluster
        }
        for (k in 1:mdp$K){
            mdp$s[i]=k              # try new cluster membership..
            if (nk[k]<ns){
                llwith[k] = mdpMarg(k,mdp) ## llwith is *with* i, llwo *w/o*
                lfyi = llwith[k]-llwo[k]
            } else 
                lfyi = mdpFy(i,k,mdp)
            lps[k] = lps[k]+lfyi
        }# k
        if (prnew>0){ # consider a new cluster..
            k=mdp$K+1
            mdp$s[i]= k  # try new cluster membership
            llwith = c(llwith, mdpMarg(k,mdp)) # incremement by 1, llwith is *with*
            lps[k] = lps[k]+llwith[k]
        }
        ps = exp(lps-max(lps))
        k = sample(1:length(ps),size=1,prob=ps)
        mdp$s[i] =k
        if (k==mdp$K+1){
            mdp$K = mdp$K+1
            out = mdpPBeta(mdp,k,T) # sample betas[k]
            mdp$betas = cbind(mdp$betas, out$betask)
            nk = c(nk,1)
            llwo = c(llwo, llwith[k])
        } else {
            nk[k] = nk[k]+1
            llwo[k] = llwith[k] ## with the new i becomes the new llwo (the next i..)
        }
    } # it
    return(mdp)
} # mdpUpd

mdpUpdB = function(mdp,spl=T){ ###** to be updated *** ###
    ## MCMC transition probs to update etas (if updateB=T)
    ## niter iterations
    mdp$betas = matrix(0,nrow=mdp$q, ncol=mdp$K) # initializing betas of right dim
    for(k in 1:mdp$K){
        out = mdpPBeta(mdp,k,spl=spl)
        mdp$betas[,k] = out$betask
    }# k
    return(mdp)
}
        
mdpPBeta = function(mdp,k,spl=F){
        ## assumes yt is correctly built
    Ck = which(mdp$s==k)
    nk = length(Ck)
    betask = NULL # default for return value betas[k] under spl=F
    if (mdp$d==1){ # p(eta*_k | yt) for patient cluster Ck
        sumy=0
        n=0
        for (r in Ck){    ## add up sumy over all pats (=rows in y and X) in cluster Ck
            idx = !is.na(yt[r,])
            if (any(idx)){ ## any observed data for this patient
                ## note: since we drop patients & prot with all NAs there should always be..
                n = n+sum(idx)
                sumy = sumy + sum(yt[r,idx])
            }
        }
        Vi = (n/mdp$sigs+1/mdp$Sbeta)
        V=1/Vi; K=sqrt(V)
        m = V*(sumy/mdp$sigs+mdp$mbeta/mdp$Sbeta)
        if (spl)
            betask = rnorm(1,m=m,sd=sqrt(V))
        else
            betask = m
    } else { # for protein cluster Ck
        ## build the posterior moment sequentially adding prot's
        Vi = matrix(0,nrow=mdp$q, ncol=mdp$q) # initialize
        Vim = rep(0,mdp$q)
        for (cl in Ck){     ## add terms for each protein in cluster Ck
            idx = !is.na(yt[,cl])  # which are observed
            if (any(idx)){  ## any observed data for this protein?
                ## note: since we drop patients & prot with all NAs there should always be..
                Vi = Vi + t(X[idx,]) %*% X[idx,]/mdp$sigs
                Vim = Vim + 1/mdp$sigs*t(X[idx,])%*%yt[idx,cl]
            }
        }
        ## debugging ##
        ## V = solve(Vi);        m = V %*% Vim
        ## cat("V*Vim=",m[1:6],"\n yt[1:3,Ck]="); print(yt[1:3,Ck[1:min(10,nk)]])
        ## now add the prior shrinkage..
        Vi = Vi + mdp$Sibeta
        Vim = Vim + mdp$Simbeta   ## Vim "V-inv*m"
        V = solve(Vi)
        K = t(chol(V))
        m = V %*% Vim
        if (spl)
            betask = m+K %*% rnorm(mdp$q)
        else
            betask = m
    }
    return(list(m=m,V=V,K=K,betask=betask))
}

mdpBetaHat = function(Ck,Vi=NULL, Vim=NULL)
{ # returns mle for a common beta for proteins in Ck, using *raw data* y
    ## only used for debugging and initialzation
    p = ncol(X)
    if (is.null(Vi) | is.null(Vim)){ # could use Vi & Vim fr incemental updates
        Vi = matrix(0,nrow=p, ncol=p) # initialize
        Vim = rep(0,p)
    }
    for (cl in Ck){     ## add terms for each protein in cluster Ck
        idx = !is.na(y[,cl])  # which are observed
        if (any(idx)){  ## any observed data for this protein?
            ## note: since we drop patients & prot with all NAs there should always be..
            Vi = Vi + t(X[idx,]) %*% X[idx,]
            Vim = Vim + t(X[idx,])%*%y[idx,cl]
        }
    }
    V = solve(Vi)
    K = t(chol(V))
    m = V %*% Vim
    return(list(m=m, Vi=Vi, Vim=Vim, V=V))
    ## returning Vi and Vim allows easy updating by adding/substracting proteins
}
    

mdpMarg = function(k,mdp)
{
    ## marginal likelih for y*[k]=y[s==k]
    ## m=mbeta, V=sig2*I + Xk' S Xk
    ## d=1 for patients; d=2 for proteins; d=0 for data records
    ##
    ## 1. posterior p(betas[k] | yk,Xk)=N(m,V), in prep for #2 below..
    out = mdpPBeta(mdp,k) ## p(betas[k] | dta)
    m = out$m
    V = out$V
    K = out$K
    Ck = which(mdp$s==k)
    sd = sqrt(mdp$sigs)
    if (mdp$d==1){ # patients
        ## 2. use candidate formula to evaluate marginal
        lfk = 0
        for(r in Ck){ ## add up log likelihood over all patients in Ck
            idx = !is.na(yt[r,])
            if (any(idx))
                lfk = lfk + sum( dnorm(yt[r,idx], m=m, sd=sd, log=T) )
        } # for r
        lfk = lfk - 0.5*log(2*pi) - log(K)
    } else { # proteins, d=2
        lfk = 0
        mk = X %*% m   ## mean over all patients for proteins in Ck
        for(cl in Ck){ ## add up log likelihood over all prots in Ck
            idx = !is.na(yt[,cl]) # patients with data on prot cl
            ## 2. use candidate formula to evaluate marginal
            lfk = lfk + sum( dnorm(yt[idx,cl], m=mk[idx], sd=sd,log=T) )
        }
        lfk = lfk - sum(log(diag(K))) - mdp$q/2*log(2*pi)
    }
    return(lfk)
}

mdpFy = function(i,k,mdp)
{
    ## likelih for y[i] | s==k
    sd = sqrt(mdp$sigs)
    if (is.na(mdp$betas[k][1])){ ## not updated - note, betas[k][1] is indep of q
        out = mdpPBeta(mdp,k,spl=T)
        mdp$betas[,k] = out$betask
    }
    lfk = 0        # default if no data for this patient or protein..
    if (mdp$d==1){ # likelihood for patient i in cluster k
        idx = !is.na(yt[i,])
        if (any(idx))
            lfk = sum( dnorm(yt[i,idx], m=mdp$betas[,k], sd=sd, log=T) )
    } else { # for protein i in cluster k, d=2
        mk = X %*% mdp$betas[,k]   ## mean over all patients for any protein in cluster k
        idx = !is.na(yt[,i]) # patients with data on prot i
        if (any(idx))
            lfk =sum( dnorm(yt[idx,i], m=mk[idx], sd=sd,log=T) )
    }
    return(lfk)
}

updateSigs = function(mdp){
    ## update mdp$sig2, residual var
    ## assume yt is already just residuals, corrected for
    ##     protein *and* patient effects (done with mdpOffset(.) before)
    S = sum(yt*yt, na.rm=T)
    N = sum(!is.na(yt))
    a = (prior$asig+N)/2
    b = (prior$bsig+S)/2
    w = rgamma(1,a,rate=b)
    sigs = 1/w
    ## save for debugging & summaries
    SSM <<- S/N

    return(sigs)
}


#######################################################
## plots and debugging

pltInit = function()
{ ## in preparation for plt(), read in all simulation summaries
    mcmc <<-  read.csv("chain.txt",header=T)
    Ey <<-  as.matrix( read.csv("Ey.txt", header=F) )
    np=ncol(Ey)
    my <<-  matrix(apply(Ey,1,mean),ncol=2) ## avg expression by condition (not used..)
    My = apply(my,1,mean)                   ## overall avg profile (never used..)
    Ey <<-  array(Ey,dim=c(16,2,np))    ## for easier access below
    ## dimensions are: age, case, protein
    
    Ey2 <<-  as.matrix( read.csv("Ey2.txt", header=F) ) 
    my2 <<-  matrix(apply(Ey2,1,mean),ncol=2)
    Ey2 <<-  array(Ey2,dim=c(16,2,np))
    
    Eyp <<-  as.matrix(read.csv("Eyp.txt", header=F) ) 
    mp <<-  matrix(apply(Eyp,1,mean),ncol=2) ## avg by condition
    Eyp <<-  array(Eyp,dim=c(16,2,np))       # Eyp as 3-d array
    mp3  <- array(mp,dim=c(16,2,np))         ## avg by condition repeated as needed..
    Eyp0 <<- Eyp-mp3        ## Eyp corrected by avg per patient (non likelihood identifiable..)
        
    Eyp2 <<- as.matrix(read.csv("Eyp2.txt", header=F) ) 
    mp2 <<-  matrix(apply(Eyp2,1,mean),ncol=2)
    Eyp2 <<-  array(Eyp2,dim=c(16,2,np))

    yy <<-  array(y,dim=c(16,2,np))   # in same 3-d array for easier plotting blow..
}

plt = function(fit=T, dta=F, prot=F, idx=NULL, lw=0.5, pltm=F, case=F,ctr=T, dtatype="p")
{
    ## plots fitted profiles, indicators select different plots (see below for examples)
    ## fit:  plot fitted mean functions for the selected proteins
    ## dta:  add dots for the observed data
    ## prot: plot fitted mean functions using protein-effects only (w/o patient effects)
    ## pltm: add average mean function, averaging across all proteins
    ## case: show selected plot for cases
    ## ctr:  .. for controls
    ##       (pltm can not be restricted to cases or controls only)
    
    ## fitted mean functions for each patient: 
    par( mar=c(4.1,4.1,0.5,0.5), mgp=c(1.5,0.5,0))
    
    np=ncol(Ey)
    npat = nrow(Ey)
    pltmch = ifelse(pltm,"l","n")
    
    if (is.null(idx))
        idx=sample(1:np,10,replace=F)
    if (fit | dta){
      ## fitted lines
      I=length(idx)
      matplot(ages, my ,type=pltmch,bty="l",col=1:2,lwd=3,ylim=range(Ey,na.rm=T), ylab="Y")
      if (fit){
        if (ctr) matlines(ages[,1], Ey[,1,idx], col=1:I, lwd=lw, lty=1, type="l")
        if (case) matlines(ages[,2], Ey[,2,idx], col=1:I, lwd=lw, lty=2, type="l")
      } # fit
      if (dta){
      ## add data
          if (ctr) matlines(ages[,1], yy[,1,idx], col=1:I, lwd=0.53, lty=1,
                            type=dtatype,pch=20)
          if (case) matlines(ages[,2], yy[,2,idx], col=1:I, lwd=0.53, lty=2,
                             type=dtatype,pch=1)
        } #dtad
    } # fit | dta
    if (prot){
      ## proteins w/o pat effects
      matplot(ages, mp ,type=pltmch,bty="l",col=1:2,lwd=3,ylim=range(Eyp,na.rm=T))
      if (ctr) matlines(ages[,1], Eyp0[,1,idx], col="grey", lwd=1, lty=1, type="l")
      if (case) matlines(ages[,2], Eyp0[,2,idx], col="pink", lwd=1, lty=2, type="l")
    } # prot
}

maxDiff = function()
{ ## reports proteins with max diff across ages
    Tm = nrow(ages)
    df = abs( (Eyp0[Tm,1,]-Eyp0[1,1,])-(Eyp0[Tm,2,]-Eyp0[1,2,]))
    of=order(-df)
    ## same with raw dta
    dfy = abs( (yy[Tm,1,]-yy[1,1,])-(yy[Tm,2,]-yy[1,2,]))
    oy=order(-dfy)
    ## with protein-specific regression
    sigs = mean(mcmc[,3])  # posterior mean for sig2
    yhat = matrix(0,nrow=32,ncol=np) # initilaize
    Vi = t(X)%*% X         # same for all proteins
    V = solve(Vi)
    beta = rep(0,p)
    for (j in 1:np){
        Vim = t(X)%*%yt[,j]
        beta = V %*% Vim
        yhat[,j] = X %*% beta
    }
    yyhat = array(yhat,dim=c(16,2,np))   # in same 3-d array for easier plotting blow..
    dfyhat = abs( (yyhat[Tm,1,]-yyhat[1,1,])-(yyhat[Tm,2,]-yyhat[1,2,]))
    oyhat=order(-dfyhat)

    nm = colnames(y)
    cat("Proteins with largest estimated effect (L2):\n  ", nm[of[1:10]],"\n")
    cat("  ", nm[of[11:20]],"\n")
    cat("Proteins with largest empirical effect (L2):\n  ", nm[oy[1:10]],"\n")
    cat("Proteins with largest mle effect (L2):\n  ", nm[oyhat[1:10]],"\n")
    return(of)
}





paper = function()
{ # calls all the plots for the paper
    pltInit()
    
    of=maxDiff()

    ps("201") # for some randomly selected patients
    plt(T,T,F,201,case=T,ctr=T)
    devoff()
    ps("1201")
    plt(T,T,F,1201,case=T,ctr=T)
    devoff()
    ps("2201")
    plt(T,T,F,2201,case=T,ctr=T)
    devoff()
    ps("15")
    plt(T,T,F,15,case=T,ctr=T)
    devoff()
    ps("115")
    plt(T,T,F,115,case=T,ctr=T)
    devoff()
    ps("protP") # plot some proteins, w/o patient effects
    plt(F,F,T,of[1:200],case=T,ctr=T,pltm=F)
    devoff()

    ps("dta")
    plt(F,T,F,of[1:20],case=T,ctr=T,dtatype="l")
    devoff()
}


main(1000)
