# SEP-BNP functions 

#######################################################
## model
##
## subjects j=1..n, cluster indicators S[j] 
##                  clusters C[k], k=1..K
## OTUs,    b=1..B, cluster indicators M[b,k]
##                  clusters, D[k,l], l=1...L (nested in C[k])
##          Note: in the notes the OTU index is "i" (=b :-)
## cluster-specific pars {th[l]=(mu[l], sig2[l]); l=1..L[k]
##                      (common atoms across k)
## data:    y[bj], B x n matrix, row=OTU, col=subj
##
## p(S,M,pi,w,mu,sig2) = prod_{j,b} p(y[bj] | th[M[i,k=S[j]]] ) *
##                                  likel
##    * prod_j pi[S_j] prod_k{ prod_i w[k, M[i,k]] } *
##       subj clusters    OTU clusters (nested wi)
##    * p(pi) * prod_k {p(w[k,.] * prod_l p(th[l])
##
##
###################################################

###################################################
## initializing pars for MC
## 

init.pi <- function()
{## pi[k], cluster probs for S[j]; generate from the prior
  vk=c(rbeta(K-1,1,alpha),1) # fractions of stick-breaking
  pi=c(vk[1],sapply(2:K,function(i) vk[i]*prod(1-vk[1:(i-1)])))
  return(pi)
} 

init.w <- function()
{## weights for G[k], k=1..K, generate from the prior
  vkl=matrix(0,K,L)          # fractions of stick-breaking
  for(k in 1:K){
    vkl[k,]=c(rbeta(L-1,1,beta),1)
  }
  w=matrix(0,K,L)
  for(k in 1:K){
    w[k,]=c(vkl[k,1],sapply(2:L,function(i) vkl[k,i]*prod(1-vkl[k,1:(i-1)])))
  }
  return(w)
}

init.Sj  <- function(K0=K)
{#initialize Sj with hclust
  ty=data.frame(t(y))
  ds=dist(ty,method='euclidean')
  hc <- hclust(ds, method = "complete" )  ## clustering subjects
  sub_grp <- cutree(hc, k = K0)
  table(sub_grp)
  ty %>% mutate(cluster = sub_grp) -> ty
  Sj=ty$cluster
  Sj0 <<- Sj  # remember intial values
  return(Sj)
}

init.muSig <- function()
{### initialize mu,sig with cluster-specific moments using the
  ### clusters under k-means (for OTU's!! Not for subjects)
  ydf <- data.frame(y)
  km  <- kmeans(ydf, centers=L, iter.max = 100, nstart = 10)  ## use k-means for initial clustering of OTU's
  ydf %>% mutate(cluster = km$cluster) -> ydf
  mu=rep(0,L)
  sig2=rep(0,L)
  for(l in 1:L){
    mu[l]= mean(as.matrix( 
      ydf[ydf$cluster == l, ] %>% dplyr::select(-cluster)
      ))
    sig2[l]= var(as.vector( 
      as.matrix(ydf[ydf$cluster == l, ]%>% dplyr::select(-cluster))
      ))
  }
  ## sig2=rep(var(as.vector(as.matrix(ydf[ydf$cluster == which.max(table(ydf$cluster)), ]%>% select(-cluster)))),L) 
  return(list(mu=mu,sig2=sig2))
}


init.Mki <- function(mu,sig2,Sj,w)
{  #same as update M_nik, except that no current Mki is used
  ## and no marginalization w.r.t mu for nl[l]=1
  nk=sapply(1:K,function(k) sum(Sj==k))
  mki = matrix(0,nrow=K,ncol=B)
  for(k in 1:K){
    zk=y[,Sj==k]
    if(nk[k]==0){
      mki[k,]= sample(1:L,B,replace=T,prob=w[k,])
    }else{
      ## mki[k,c(1)]=1
      for(b in 1:B){
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
        }# for l
        log_prob=log_prob-max(log_prob)
        mki[k,b]=sample(1:L,1,prob = exp(log_prob))
      }# b
    }# else
  }# k
  return(mki)
}

#######################################################
## evaluating log likelihood & posterior
## (for debugging, and for MH moves)
##

ll_fun <- function(mu,sig2,Sj,mki,w=NULL,pi=NULL,post=F){
  ## loglikelihood check functions
  ## if post, then return log posterior (needs additional arg's)
  ll=0; lpr=0
  for(j in 1:n){
    k=Sj[j]
    mkik = mki[k,]
    ll = ll+sum(dnorm(y[,j], mu[mkik], sd=sqrt(sig2[mkik]), log=T))
  }
  if (post)
    lpr = lprior(mu,sig2,Sj,mki,w,pi,ll)
  return(ll+lpr)
}

## for debugging only -- keepin track of each factor in log posterior
## in the two most recent (old, new) posterior evaluations 
lpot = c("lik","mu","sig","Sj","mki","pi","w")
lpo = matrix(0,nrow=2,ncol=length(lpot)) ## to save most recent 2 eval of log post
colnames(lpo)=lpot
rownames(lpo)=c("old","new")

lprior  <- function(mu,sig2,Sj,mki,w,pi,ll)
{ # log prior; argument "ll" only for debugging..
  lpmu = sum( dnorm(mu, m=mu0, sd=sig0, log=T) )
  lpsig2 = sum( dinvgamma(sig2, a0, b0, log=T) )
  lpSj = sum( log(pi[Sj]) )
  lpmki = 0
  for(k in 1:K)
    lpmki = lpmki + sum( log(w[k,mki[k,]]) )
  vk = pi
  for(k in 2:K)
    vk[k] = pi[k]/prod(1-vk[1:(k-1)])
  lpi = sum( dbeta(vk[1:(K-1)], 1, alpha,log=T) )
  vkl = w
  lw=0
  for(k in 1:K){
    for(l in 2:L)
      vkl[k,l] = w[k,l]/prod(1-vkl[k,1:(l-1)])
    lw = lw + sum( dbeta(vkl[k,-L], 1, beta, log=T) )
  }# k
  ## for debugging only - save all factors of two most recent
  ## log *post* evaluations (including ll, which is passed as arg)
  lpo[1,] <<- lpo[2,] 
  lpo[2,] <<- c(ll, lpmu, lpsig2, lpSj,lpmki,lpi,lw)
  
  lpr = lpmu+lpsig2+lpSj+lpmki+lpi+lw
  return(lpr)
}



check.ll <- function(txt,lastll,mu,sig2,Sj,mki,w=NULL, pi=NULL, post=F)
{## for debugging - check if log posterior drops by
  ## more than 500 -- should never happen! (this is a bit *very* conservative
  ll  <- ll_fun(mu,sig2,Sj,mki,w,pi,post)
  if (ll < lastll-500) {## bug!
    cat("after ",
        txt,"\t ll=", format(ll), "\t last-ll=",format(lastll),"\n")
    print(lpo)
    browser()
  }
  return(ll)
}




#######################################################
## MCMC transition prob's 

Update_mu <- function(yl,nl,sig2_l){
  ## =3. in the notes
  if (nl==0) # empty OTU cluster - generate from prior
    ##        mu = rtruncnorm(1,0,quantile(as.matrix(y),0.9),mu0,sig0)
    mu = rnorm(1,mu0,sd=sig0)
  else{
    v=1/(1/sig0^2 + nl/sig2_l)
    m=v*(mu0/sig0^2+sum(yl)/sig2_l)
    ## mu=rtruncnorm(1,0,Inf,m,sd=sqrt(v))
    mu=rnorm(1,m,sd=sqrt(v))
  }
  return(mu)
}


Update_sig2 <- function(yl,nl,mu_l){
  ## =4. in the notes
  if (nl==0) # empty OTU cluster - generate from prior
    sig2=rinvgamma(1,a0,b0) 
  else
    sig2=rinvgamma(1,a0+nl/2,b0+sum((yl-mu_l)^2)/2)
  return(sig2)
}


update.pi <- function(Sj)
{  # update vk and pi where pi[k] = v[k]*prod{h<k} (1-v[h])
  ## =6. in the notes
  nk=sapply(1:K,function(k) sum(Sj==k))
  vk=c(sapply(1:(K-1),function(i)
    rbeta(1,1+nk[i],alpha+sum(nk[(i+1):K]))),1)
  pi=c(vk[1],sapply(2:K,function(i) vk[i]*prod(1-vk[1:(i-1)])))
  return(pi)
}

update.w <- function(mki)
{ #update w_k
  ## =5. in the notes
  ## w[k,l] = v[k,l]*prod{h<l} (1-v[k,h])
  vkl = matrix(0,nrow=K,ncol=L)
  w = 0*vkl
  for(k in 1:K){
    nkl=sapply(1:L,function(l) sum(mki[k,]==l))
    vkl[k,]=c(sapply(1:(L-1),function(i)
      rbeta(1,1+nkl[i],beta+sum(nkl[(i+1):L]))),1)  
    w[k,]=c(vkl[k,1],sapply(2:L,function(i) vkl[k,i]*prod(1-vkl[k,1:(i-1)])))
  }
  return(w)
}


update.SjMarg <- function(mu,sig2,Sj,pi,w)
{ #update Sj using MH trans prob
  ## TQARGET: Cond prob p(Sj=k | ...) requires eval with all subj's
  ## p(Sj=k | ...) \propto                                  eq (2) in the notes
  ## prod_k"
  ##   prod_i
  #      { sum_l
  ##       { prod_{j" in C[k"]} N(y[ij"] | th[l])
  ##       }
  ##     }
  ## note: either (A) evaluate for all k" (and can simplify a bit) (Gibbs) (as in eq (2))
  ##          or  (B) evaluate (*) for propoed k and current kold (MH) (eq (3) in notes)
  ##       We will do (B) (but could just as well do (A))
  ## MH PROPOSAL:
  ##   q(Sj=k) \propto pi[k] prod_i {sum_l w[kl] N(y[ik] | th[l]) != (*)
  ## PLAN: following (B) we
  ## (i) generate proposal 'k' using q (evaluating q[k] for k=1...K
  ## (ii) MH acc prob, evaluating (*) then only for proposal 'k' and
  ##      current 'kold'
  lq = rep(0,K) # initialize
  lpyk = rep(0,4) ## just save p(y*[k] | Sj=k..) for debugging..
  ## lpyk[1]=p(y[*,j] | Sj=current kold);
  ##     [2]=              new     k)
  ##     [3]=
  for(j in 1:n){
    ## (i) MH proposal using just the mixture model Gk for y[j]
    ##    (which is distinct from the mixture of mv normals (*))
    kold = Sj[j]           ## remove from current cluster
    for (k in 1:K)
      lq[k] = log(pi[k])+
        sum( sapply(1:B,
                    function(i) dGk(y[i,j],mu=mu,sig2=sig2,wk=w[k,])))
    ## use this as proposal prob's
    lq = lq - max(lq)
    k=sample(1:K,1,prob=exp(lq))
    if (k==kold) next   ## nothing changes..
    ## (ii) Acc ratio; now we evaluate (*) for clusters (k,kold)
    ##      once with Sj=k, and once for Sj=kold
    lpyk[1] = dGkmv(kold,mu,sig2,Sj,w)  # current, using direct comp
    lpyk[2] = dGkmv(k   ,mu,sig2,Sj,w)
    dpyold =  lpyk[1]+lpyk[2]
    Sj[j]=k
    lpyk[3] = dGkmv(kold,mu,sig2,Sj,w)  # proposal, using direct comp
    lpyk[4] = dGkmv(k   ,mu,sig2,Sj,w)
    dpy       =  lpyk[3]+lpyk[4]      
    A = (dpy+log(pi[k]) - (dpyold+log(pi[kold])))+(lq[kold]-lq[k])
    if (log(runif(1)) > A){ # reject
      Sj[j]=kold
    }
  }# j
  return(Sj)
}

## ********************************* aux functions for update.Sj

dGk <- function(y,mu,sig2,wk)
{# evaluates *log* G_k(y), marginalizing mki, conditioning on mu,sig,w
  ## careful! THis is the marginal for p(y[j] | Sj=k,mu,sig2,pi)
  ## but when margalizing wrt M[ik] the y[ij], j \in {j": S[j"]=k} 
  # become dependent!
  pyl = dnorm(y, mu, sd=sqrt(sig2), log=F)
  py = sum(wk*pyl)
  lpy = log(py)
  return(lpy)
}

dGkmv <- function(k,mu,sig2,Sj,w)
{# evaluates *log* p(y*[k,] | mu,sig2 )
  ## marginalizing mki, conditioning on mu,sig,w
  ## see the comment about marginal vs. joint above
  jdx = which(Sj==k)
  nk = length(jdx)
  yk = as.matrix( y[,jdx] )
  lpybl = rep(0,L)
  mx = rep(0,B)
  lpy = 0
  for(b in 1:B){
    for(l in 1:L){
      lpybl[l] = sum( dnorm(yk[b,], mu[l], sd=sqrt(sig2[l]), log=T) )
    }# l
    mx[b] = max(lpybl)
    pyb   = sum(w[k,]*exp(lpybl-mx[b]))  ## * e^mx
    lpy   = lpy+ log(pyb)
  }# b
  lpy = lpy + sum(mx)
  return(lpy)
}

## *******************************************************



update.Mki <- function(mu,sig2,Sj,w,mki)
{  #update M_nik
  ## = item 2. in notes
  nk=sapply(1:K,function(k) sum(Sj==k))
  nl = rep(0,L)
  for(k in 1:K)
    nl=nl+nk[k]*sapply(1:L,function(l) sum(mki[k,]==l))
  oldmki = mki # remember for counting # changes later..
  
  for(k in 1:K){
    if(nk[k]==0){
      mki[k,]= sample(1:L,B,replace=T,prob=w[k,])
    }else{
      yk=as.matrix( y[,Sj==k] ) ## as.matrix needed for nk=1
      for(b in 1:B){
        oldl = mki[k,b]    # remove current allocation
        nl[oldl] = nl[oldl]-nk[k]
        log_prob=rep(-Inf,L)
        for(l in 1:L){ # loop over all non-zero OTU-clusters
          if (nl[l]==0){ # use marginal, wrt mu[l]
            ybk = yk[b,]
            sum_loglk=lmy.Mki(ybk, sig2[l])
          } else {
            m=mu[l]; s= sqrt(sig2[l])
            sum_loglk=sum( dnorm(yk[b,],m,s,log=T) )
            if(sum_loglk==-Inf)  sum_loglk=-10^5
          }
          log_prob[l]=log(w[k,l])+sum_loglk
        }# for l
        log_prob=log_prob-max(log_prob)
        l = sample(1:L,1,prob = exp(log_prob))
        mki[k,b]= l
        nl[l] = nl[l]+nk[k]
        if (nl[l]==nk[k]){ ## opening new cluster
          mu[l]=Update_mu(yl=ybk,nl=nl[l], sig2_l=sig2[l])
        }# new cluster
      }# b
    }# else, nk[k]>0
  }# k
  nmki = sum( oldmki[nk>0,] != mki[nk>0, ] )  # n (non-trivial) changes
  Nmki = B*sum(nk>0)                          # n transitions
  return(list(mki=mki,nmki=nmki,Nmki=Nmki,mu=mu))
}

## ********************************* aux fct for update.mki
lmy.Mki  <- function(ybk,sig2l)
{## marginal p(y*[b,k] | Mki[k,b]=l, nl[l]=0); y*[b,k]=ybk, nl=length(ybk)
  ## use candiate formula
  ## p(y*[b,k] |...) = p(y*[b,k] | mu[l] ...)*p(mu[l])/p(mu[l] | y*[b,k]...)
  ## using posterior mean mu[l] = E(mu[l] | y*[b,k]...) as in update_mu
  ## mu[l] = v*(mu0/sig0^2 + sum(y*[b,k]/sig2_l))
  ##         1/v = 1/sig0^2 + nl/sig2_l,  (nl = nk[k])
  ## p(mu[l] | ...) = (2*pi*v)^{-1/2}
  ## p(y*[b,k] | mu) = prod N(y[b,j] | mu,sig2_l), Sj=k
  ## p(mu[l]) = N(mu0,sig0^2)
  nkk = length(ybk)
  v = 1/(1/sig0^2 + nkk/sig2l)
  m = v*(mu0/sig0^2  + sum(ybk/sig2l))
  lym = sum( dnorm(ybk, m, sd=sqrt(sig2l), log=T) )
  lm = dnorm(m, m=mu0, sd=sig0, log=T)
  lmy = -1/2*log(2*3.141593*v)
  ly = lym+lm-lmy
  return(ly)
}


#######################################################
## mcmc

writeMCMC <- function(iter,pi,Sj,w,mki,mu,sig2,ll,pmki=rep(0,L),app=T)
{ # open files to save MCMC sims
  ## pmki = # accepted moves for each of the L OTU clusters
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
}

niter=500; iter=0
mcmc  <- function(niter=1000, pi, w, Sj, mki, mu, sig2, niter0=500)
{#### MCMC
  ## estimate Sj after niter0 iterations, and stop updating Sj
  last_ll=check.ll("update_w", lastll= -1e10, mu,sig2,Sj,mki,w,pi,post=T)-100000
  writeMCMC(0,pi=pi,Sj=Sj,w=w,mki=mki,mu=mu,sig2=sig2,ll=last_ll,app=F,pmki=0)
  nmki = 0  # keep track of n of moves of mki
  Nmki = 0  # n of transitions (=iter * K+ * L)
  SjMC = NULL
  add  = F
  for(iter in 1:niter){
    pi = update.pi(Sj=Sj)
    last_ll=check.ll("update_pi",last_ll, mu,sig2,Sj,mki,w,pi,post=T)
    w = update.w(mki=mki)
    last_ll=check.ll("update_w",last_ll, mu,sig2,Sj,mki,w,pi,post=T)
    if (iter<niter0){
      Sj = update.SjMarg(mu=mu,sig2=sig2,Sj=Sj,pi=pi,w=w) # marg wrt Mki[nk==0,]
      ## last_ll=check.ll("update_Sj",last_ll, mu,sig2,Sj,mki) mki is invalid now!
    }
    out =
      update.Mki(mu=mu,sig2=sig2,Sj=Sj,w=w,mki=mki) # marg wrt mu[nkl==0]
    mki = out$mki
    mu = out$mu
    nmki = nmki+out$nmki
    Nmki = Nmki+out$Nmki
    last_ll=check.ll("update_Mki",last_ll, mu,sig2,Sj,mki,w,pi,post=T)
    
    ## update mul and sigl
    mki_expand=t(mki[Sj,]) # B x n matrix of M[b, Sj]
    checkMki(mki,mki_expand,Sj)
    for(l in 1:L){
      yl= y[mki_expand==l]
      nl=sum(mki_expand==l)
      mu[l]=Update_mu(yl=yl,nl=nl,sig2_l=sig2[l])
      last_ll=check.ll(paste("update_mu",l),last_ll, 
                       mu,sig2,Sj,mki,w,pi,post=T)
      sig2[l]=Update_sig2(yl=yl,nl=nl,mu_l=mu[l])
      last_ll=check.ll(paste("update_Sig2",l),
                       last_ll, mu,sig2,Sj,mki,w,pi,post=T)
    } #
    txt <- paste("it ",format(iter))
    last_ll=check.ll(txt, last_ll, mu=mu,sig2=sig2,Sj=Sj,mki=mki,w,pi,post=T)
    if (iter %% 10 == 0){ # save every 10th iter
      pmki = nmki/Nmki
      writeMCMC(iter, pi, Sj, w, mki, mu, sig2=sig2, ll=last_ll, pmki=pmki)
      plt.Gk(w, mu, sig2, add, Sj)
      add=F
      print(iter)
    }
    if (iter < niter0)
      SjMC = rbind(SjMC, Sj) # for next step...
    if (iter==niter0){ # estimate Sj
      Sj = salso(SjMC)
      cat("\n Fixing Sj using 'salso'.\n")
      ## relabel subj clusters by decreasing size
      nk=sapply(1:K,function(k) sum(Sj==k))
      Sidx = order(nk, decreasing=T)
      Srk = rank(-(nk), ties="first") ## runif to avoid ties
      Sj = Srk[Sj]
      ## w,pi,mki are totally out of synch now & need re-init:
      ## impute new pi, mki
      pi = update.pi(Sj=Sj)
      out =
        update.Mki(mu=mu,sig2=sig2,Sj=Sj,w=w,mki=mki) # marg wrt mu[nkl==0]
      mki = out$mki
      mu = out$mu
      w = update.w(mki=mki)
      last_ll = -10e10  # reset last ll
    }
  }# iter
  return(Sj)
}

## SjMC = as.matrix(read.table("Sj.txt"))
salso <- function(C){
  ## finds point estimate for random partition based on MC sample C
  ## (each row of C is a random partition, represented by cluster memberhsip
  ##         indicators)
  ## C = (M x J), with C[m,] = cluster membership indicators for
  ##                        partition of J units
  ## Find cstar to minimize posterior expected loss
  ## mstar = arg min_m \sum_{m2} d(C[m,], C[m2,]) \approx
  ##                   E{ d(C[m,], c]) | y }
  ## cstar = C[mstar,]
  
  M=nrow(C)
  J=ncol(C)
  dmm = matrix(0,nrow=M,ncol=M) #d(Sj[m], Sj[m'])
  ##
  A  <- NULL # A[m,] = adjacency matrix for m-th sample
  for (m in 1:(M-1)){ # adjacency matrix for m-th sample
    Am = ifelse(dist(C[m,])>0,0,1)
    A <- rbind(A,Am)
    if (m==1) next
    for(m2 in 1:(m-1))
      dmm[m,m2] = sum(abs( A[m,]-A[m2,] )) # hamming distance d(Am,Am2)
  }# m
  dmm = dmm+t(dmm) # fill in upper triangular
  Lm = apply(dmm,1,sum)
  ms = which.min(Lm)[1]
  return(C[ms,])
}

## just for debugging...
checkMki  <- function(mki, mki_expand, Sj)
{ 
  for(j in 1:n){
    for(b in 1:B){
      if (mki_expand[b,j] != mki[Sj[j],b] ){
        browser()
      }
    }
  }
}

#######################################################
## main driver
ex <- function(niter=2000, niter0=250, niter1=500){
  set.seed(1992)
  pi = init.pi()
  w = init.w()
  Sj = init.Sj(3)
  out = init.muSig()
  mu = out$mu
  sig2 = out$sig2
  mki = init.Mki(mu,sig2,Sj,w)
  
  Sj = mcmc(niter, pi=pi,w=w,Sj=Sj,mki=mki,mu=mu,sig2=sig2, niter0=niter0)
  klist = fig4(niter,niter0,niter1)
  # rcR(1,2)
  fig5(klist, niter,niter0,niter1)
} 

#######################################################
## make figures,  also called at the end of ex(),
## but could be called separately

niter = 1000; niter0 = 250; niter1 = 750
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
  idx <- order(s, decreasing=T)
  kdx  <- which( nk>0 )
  sk = matrix(0, nrow=K, ncol=B)
  for (k in klist){ # non-empty clusters; clusters k=1,2,3
    sk[k,] = apply(Y[,Sj==k], 1, sum)
  }
  Sk = apply(sk, 1, sum)
  sk = sk/Sk
  csk = apply(sk[,idx], 1, cumsum) 
  # cum sum for each k, OTU's ordered by total frequ
  lwd = nk/max(nk)*3
  matplot(1:B, csk[,klist], type="l", xlab="OTU", ylab="CUMSUM",
          lwd=lwd[klist], bty="l")
  legend(60, 0.6, col=klist, lty=klist, legend=klist, bty="n")
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

# Reorder rows and columns (observations) of a dissimilarity matrix intra groups 
# and possibly reorder also the groups (batch of observations)
reorder_dismat <-  function(dismat, groups, order.groups=NULL){
  # Use correlation between variables as distance
  order.dis   = integer(0)
  J           = length(unique(groups))
  if(is.null(order.groups)){
    order.j   = 1:J
  } else {
    order.j   = order.groups
  }
  for (j in order.j){
    groups.j  = which(groups==j)
    dd        = as.dist((1-dismat[groups.j, groups.j])/2)
    hc        = hclust(dd)
    order.dis = c(order.dis, hc$order+length(order.dis))
  }
  dismat      = dismat[order.dis, order.dis]
  dismat      = dismat[nrow(dismat):1,]
}

## Function to plot the heatmap of the posterior probabilities of co-clustering
## of obs assigned to vertices
Plot_heat <- function(dissimlar_stable = dissimlar_stable,
                      I          = B){
  dismat      = round(dissimlar_stable, 2)
  dismat      = reorder_dismat(dismat,groups=rep(1,I))
  plot_dismat = reshape2::melt(dismat)
  ggplot(data=plot_dismat, aes(x=factor(Var1), y=factor(Var2), fill=value)) + 
    geom_tile()+ theme_bw()+ 
    scale_y_discrete(breaks = floor(seq(1, I, length.out = 9)), 
                     labels = floor(seq(1, I, length.out = 9))) +
    scale_x_discrete(breaks = floor(seq(1, I, length.out = 9)), 
                     labels = floor(seq(1, I, length.out = 9))) +
    xlab("OTU") + ylab("OTU") +
    scale_fill_gradientn(colours = c("white", "yellow", "red"), 
                         values = rescale(c(0, 0.5, 1)), 
                         space = "Lab", name="") +
    theme(legend.position = "right", text = element_text(size=20))
}
