# SEP-BNP functions 

# SEP-RPM

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

# Assign K = 10 as default if not assigned
if(!exists("K")){K<-10}
if(!exists("L")){L<-15}

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
    # print(lpo)
    # browser()
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
## MCMC
writeMCMC <- function(iter,pi,Sj,w,mki,mu,sig2,ll,pmki=rep(0,L),app=T)
{ # open files to save MCMC sims
  ## pmki = # accepted moves for each of the L OTU clusters
  write(pi, "Data-and-Results/pi.txt", K, append=app)
  write(Sj, "Data-and-Results/Sj.txt", n, append=app)
  write(c(w), "Data-and-Results/w.txt", L , append=app)
  write(c(mki), "Data-and-Results/mki.txt", K*B, append=app)
  write(mu, "Data-and-Results/mu.txt", L, append=app)
  write(sig2, "Data-and-Results/sig.txt", L, append=app)
  nk=sapply(1:K,function(k) sum(Sj==k))
  write(nk, "Data-and-Results/nk.txt", K, append=app)
  nl=sapply(1:L,function(l) sum(c(mki)==l))
  write(nl, "Data-and-Results/nl.txt", L, append=app)
  write(ll, "Data-and-Results/logl.txt", 1, append=app)
  
  ## write summary of the iteration to "iter.txt"
  nk=sapply(1:K,function(k) sum(Sj==k))
  nl= sapply(1:L, function(l) sum(c(mki)==l))
  nkp=sum(nk>0)
  nlp=sum(nl>0)
  summ = c(iter, nkp, nlp, nk, nl, ll, pmki)
  write(summ, "Data-and-Results/iter.txt", length(summ), append=app)
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
  SjMC  <- read.myfile("Data-and-Results/Sj.txt",n)
  mu = read.myfile("Data-and-Results/mu.txt", L)
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
  mki2 <- read.myfile("Data-and-Results/mki.txt",K)
  summ = read.myfile("Data-and-Results/iter.txt")
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
fn <- "Data-and-Results/w.txt"; p <- K; q  <-  L
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
  w  <- read.myfile("Data-and-Results/w.txt",K,L) # w[k, iter, l]
  niter  <- dim(w)[2]
  pi <- read.myfile("Data-and-Results/pi.txt",K)
  mu  <- read.myfile("Data-and-Results/mu.txt",L)
  sig2  <- read.myfile("Data-and-Results/sig.txt",L)
  Sj  <-  read.myfile("Data-and-Results/Sj.txt",K)
  
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

#### Functions for protein regression
readDta_reg = function(file="data_protein.RData"){
  ## cleaning data and construct design matrix
  ## rownames(dta)=c()
  load(file)
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

main_reg = function(niter=100){
  ## main driver
  sigs = 1/ (prior$asig/prior$bsig) ## prior mean for residual var
  mdpEta = mdpInit_reg(2,C,q=p, mbeta=prior$meta, Sbeta=prior$Seta,
                   a=prior$aeta, b=prior$beta, sigs=sigs)
  mdpXi =  mdpInit_reg(1,R,q=1, mbeta=prior$mxi, Sbeta=prior$Sxi,
                   a=prior$axi, b=prior$bxi, sigs=sigs)
  ## initialize partitions using hclust;
  ## in mdpInit_reg they were initialized with K=1 cluster
  mdpOffset_reg(mdpEta)
  mdpXi  = mdpInitClust_reg(mdpXi, R) # all singleton clusters
  mdpOffset_reg(mdpXi)
  Kmax = if(mdpEta$b>0) round( mdpEta$a/mdpEta$b ) else 20
  mdpEta = mdpInitClust_reg(mdpEta, Kmax)
  
  mcmcInit_reg()
  for(it in 1:niter){
    mdpOffset_reg(mdpXi)  # updates global yt = y - pat effects
    mdpEta = mdpUpdate_reg(mdpEta,niter=500)
    mdpOffset_reg(mdpEta)  # updates global yt = y - prot effects
    mdpXi = mdpUpdate_reg(mdpXi)
    mdpOffset_reg(mdpXi,incr=T)  # yt = y - (pat effects + prot effects)
    mdpXi$sigs = updateSigs_reg(mdpXi)
    mdpEta$sigs = mdpXi$sigs    # need residual var in both
    if (it %% 10==0)
      mcmcUpd_reg(mdpEta,mdpXi,it)
  }# it
}

mcmcInit_reg = function(){
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

mcmcUpd_reg = function(mdpEta, mdpXi,it){
  ## updates global summaries and writes to file (if desired..)
  ## prints out summary at iteration
  
  if(F){
  ## print summaries
  cat(it, "\t K-Prot=",mdpEta$K," (",table(mdpEta$s),"), \n avg(betas)=",
      format(apply(mdpEta$betas,1,mean),digits=2),"\n")
  cat("\t K-Pat= ",mdpXi$K," (",table(mdpXi$s),") avg(betas)=",
      format(mean(mdpXi$betas),digits=2),"\t")
  cat("sig=", format(sqrt(mdpXi$sigs)),"   SSM=", format(SSM), "\n")
  }
  
  ## summaries
  nkProt = sort( table(mdpEta$s), dec=T)[1:5]
  nkPat  = sort( table(mdpXi$s),  dec=T)[1:5]
  nkProt = ifelse(is.na(nkProt),0,nkProt) # replace NA's by 0
  nkPat = ifelse(is.na(nkPat),0,nkPat)
  
  line = c(it, SSM, mdpXi$sigs, mdpEta$K, mdpXi$K, nkProt, nkPat)
  names(line) = c("it", "SSM", "sig2", "K-prot", "K-pat", paste("nk",1:5,sep=""), paste("nk",1:5,sep=""))
  chain <<- rbind(chain,line)
  ## it, SSM, sig2, K-pat, K-prot, nk-pat*, nk-prot* (nk*: 5 largest cluster sizes, padded with 0's if needed)
  
  ## updated fitted summaries
  yhat = mdpFitted_reg(mdpEta, mdpXi)
  yphat = mdpFitted_reg(mdpEta, fitProt=T)
  Ey <<- Ey+yhat
  Ey2 <<- Ey2+yhat*yhat
  Eyp <<- Eyp+yphat
  Eyp2 <<- Eyp2+yphat*yphat
  nEy <<- nEy+1
  # browser()
  if (it %% 50 == 0){
    options(digits=2)
    write.table(chain,"Data-and-Results/chain.txt",sep=",",append=F, 
                col.names=!file.exists("Data-and-Results/chain.txt"))
    write.table(mdpXi$s,"Data-and-Results/sPat.txt",sep=",",append=filesAppend, 
                col.names=!file.exists("Data-and-Results/sPat.txt"))
    write.table(mdpEta$s, #[1:250],
                "Data-and-Results/sProt.txt",sep=",",append=filesAppend,
                col.names=!file.exists("Data-and-Results/sProt.txt"))
    ## Save all, otherwise save only the first 250 proteins if it is enough
    filesAppend <<- TRUE # from now on append those three files..
    ## the files below are always overwritten
    write.table(format(Ey/nEy), "Data-and-Results/Ey.txt", quote=F,
                col.names=F, row.names=F, sep=",")
    write.table(format(Ey2/nEy), "Data-and-Results/Ey2.txt",
                quote=F, col.names=F, row.names=F,sep=",")
    write.table(format(Eyp/nEy), "Data-and-Results/Eyp.txt",
                quote=F, col.names=F, row.names=F, sep=",")
    write.table(format(Eyp2/nEy), "Data-and-Results/Eyp2.txt", 
                quote=F, col.names=F, row.names=F, sep=",")
    print(it)
  }
  save(yt,file="Data-and-Results/yt.RData")
  return(0)
}

mdpFitted_reg = function(mdpEta, mdpXi=NULL, fitProt=F){
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

mdpOffset_reg = function(mdp,incr=F){
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

mdpSwap_reg = function(mdp,k1,k2){
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

mdpInitClust_reg = function(mdp,K0){
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
  mdp = mdpUpdB_reg(mdp,spl=F)
  return(mdp)
}

mdpUpdate_reg = function(mdp, niter=mdp$M){
  mdp = mdpUpdS_reg(mdp,niter)
  mdp = mdpUpdB_reg(mdp)
  return(mdp)
}

mdpInit_reg = function(d,M,q,
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

mdpUpdS_reg = function(mdp, niter=mdp$M, ns=20){
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
      llwo[k] = mdpMarg_reg(k,mdp) # initialize llwo
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
      mdp = mdpSwap_reg(mdp,sold,mdp$K)
      llwith[mdp$K] = llwo[sold] # currently not used - but save it
      llwo[sold] = llwo[mdp$K]
      nk[sold] = nk[mdp$K]
      nk = nk[-mdp$K]            # drop last element
      mdp$K = mdp$K-1
    } else {
      if (nk[sold]<ns){      # llwo was not evaluated; get it..
        llwith[sold] = llwo[sold] # currently not used - but save it
        llwo[sold] = mdpMarg_reg(sold,mdp) # evaluate (new) llwo (w/o i)
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
        llwith[k] = mdpMarg_reg(k,mdp) ## llwith is *with* i, llwo *w/o*
        lfyi = llwith[k]-llwo[k]
      } else 
        lfyi = mdpFy_reg(i,k,mdp)
      lps[k] = lps[k]+lfyi
    }# k
    if (prnew>0){ # consider a new cluster..
      k=mdp$K+1
      mdp$s[i]= k  # try new cluster membership
      llwith = c(llwith, mdpMarg_reg(k,mdp)) # incremement by 1, llwith is *with*
      lps[k] = lps[k]+llwith[k]
    }
    ps = exp(lps-max(lps))
    k = sample(1:length(ps),size=1,prob=ps)
    mdp$s[i] =k
    if (k==mdp$K+1){
      mdp$K = mdp$K+1
      out = mdpPBeta_reg(mdp,k,T) # sample betas[k]
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

mdpUpdB_reg = function(mdp,spl=T){ ###** to be updated *** ###
  ## MCMC transition probs to update etas (if updateB=T)
  ## niter iterations
  mdp$betas = matrix(0,nrow=mdp$q, ncol=mdp$K) # initializing betas of right dim
  for(k in 1:mdp$K){
    out = mdpPBeta_reg(mdp,k,spl=spl)
    mdp$betas[,k] = out$betask
  }# k
  return(mdp)
}

mdpPBeta_reg = function(mdp,k,spl=F){
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

mdpBetaHat_reg = function(Ck,Vi=NULL, Vim=NULL)
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

mdpMarg_reg = function(k,mdp)
{
  ## marginal likelih for y*[k]=y[s==k]
  ## m=mbeta, V=sig2*I + Xk' S Xk
  ## d=1 for patients; d=2 for proteins; d=0 for data records
  ##
  ## 1. posterior p(betas[k] | yk,Xk)=N(m,V), in prep for #2 below..
  out = mdpPBeta_reg(mdp,k) ## p(betas[k] | dta)
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

mdpFy_reg = function(i,k,mdp)
{
  ## likelih for y[i] | s==k
  sd = sqrt(mdp$sigs)
  if (is.na(mdp$betas[k][1])){ ## not updated - note, betas[k][1] is indep of q
    out = mdpPBeta_reg(mdp,k,spl=T)
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

updateSigs_reg = function(mdp){
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

pltInit_reg = function()
{ ## in preparation for plt_reg(), read in all simulation summaries
  mcmc <<-  read.csv("Data-and-Results/chain.txt",header=T)
  Ey <<-  as.matrix( read.csv("Data-and-Results/Ey.txt", header=F) )
  np = ncol(Ey)
  my <<-  matrix(apply(Ey,1,mean), ncol=2) ## avg expression by condition (not used..)
  My = apply(my,1,mean)                   ## overall avg profile (never used..)
  Ey <<-  array(Ey,dim=c(16,2,np))    ## for easier access below
  ## dimensions are: age, case, protein
  
  Ey2 <<-  as.matrix( read.csv("Data-and-Results/Ey2.txt", header=F) ) 
  my2 <<-  matrix(apply(Ey2,1,mean),ncol=2)
  Ey2 <<-  array(Ey2,dim=c(16,2,np))
  
  Eyp <<-  as.matrix(read.csv("Data-and-Results/Eyp.txt", header=F) ) 
  mp <<-  matrix(apply(Eyp,1,mean),ncol=2) ## avg by condition
  Eyp <<-  array(Eyp,dim=c(16,2,np))       # Eyp as 3-d array
  mp3  <- array(mp,dim=c(16,2,np))         ## avg by condition repeated as needed..
  Eyp0 <<- Eyp-mp3        ## Eyp corrected by avg per patient (non likelihood identifiable..)
  
  Eyp2 <<- as.matrix(read.csv("Data-and-Results/Eyp2.txt", header=F) ) 
  mp2 <<-  matrix(apply(Eyp2,1,mean),ncol=2)
  Eyp2 <<-  array(Eyp2,dim=c(16,2,np))
  
  yy <<-  array(y,dim=c(16,2,np))   # in same 3-d array for easier plotting
}

plt_reg = function(fit=T, dta=F, prot=F, idx=NULL, lw=0.5, pltm=F, case=F,ctr=T, dtatype="p")
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

maxDiff_reg = function()
{ ## reports proteins with max diff across ages
  Tm = nrow(ages)
  df = abs( (Eyp0[Tm,1,]-Eyp0[1,1,])-(Eyp0[Tm,2,]-Eyp0[1,2,]))
  of=order(-df)
  ## same with raw dta
  dfy = abs( (yy[Tm,1,]-yy[1,1,])-(yy[Tm,2,]-yy[1,2,]))
  oy=order(-dfy)
  # 
  np = ncol(Ey)
  #
  ## with protein-specific regression
  sigs = mean(mcmc[,3])  # posterior mean for sig2
  yhat = matrix(0,nrow=32, ncol=np) # initialize
  Vi = t(X)%*% X         # same for all proteins
  V = solve(Vi)
  beta = rep(0,p)
  for (j in 1:np){
    Vim = t(X)%*%yt[,j]
    beta = V %*% Vim
    yhat[,j] = X %*% beta
  }
  yyhat = array(yhat,dim=c(16,2,np))   # in same 3-d array for easier plotting
  dfyhat = abs( (yyhat[Tm,1,]-yyhat[1,1,])-(yyhat[Tm,2,]-yyhat[1,2,]))
  oyhat=order(-dfyhat)
  
  nm = colnames(y)
  cat("Proteins with largest estimated effect (L2):\n  ", nm[of[1:10]],"\n")
  cat("  ", nm[of[11:20]],"\n")
  cat("Proteins with largest empirical effect (L2):\n  ", nm[oy[1:10]],"\n")
  cat("Proteins with largest mle effect (L2):\n  ", nm[oyhat[1:10]],"\n")
  return(of)
}

### plot ggplot reg
plt_reg_ggplot = function(fit=T, dta=F, prot=F, idx=NULL, 
                          lw=0.5, pltm=F, case=F,ctr=T, dtatype="p"){
  np = ncol(Ey)
  npat = nrow(Ey)
  pltmch = ifelse(pltm,"line","n")
  
  if (is.null(idx)){
    idx=sample(1:np,10,replace=F)
  }
  if (fit | dta){
    ## fitted lines
    I=length(idx)
    if(pltmch=="n"){
      df= data.frame(ages, my)
      Plot = ggplot(data=df, aes(x=X1,y=X1.1))+
        coord_cartesian(ylim=range(Ey,na.rm=T), xlim=range(ages,na.rm=T))+
        ylab("Y") +xlab("ages")
    }
    if (fit){
      if (ctr){
        df = data.frame(ages = ages[,1], Ey = Ey[,1,idx])
        df = melt(df,id="ages")
        # Create the plot
        Plot <-Plot + 
          geom_line(data=df, aes(x = ages, y = value, group = variable, 
                    color = factor(variable)), size = lw, linetype = 1)+
          theme(legend.position="none")
      }
      if (case){
        df = data.frame(ages = ages[,2], Ey = Ey[,2,idx])
        df = melt(df,id="ages")
        # Create the plot
        Plot <-Plot + 
          geom_line(data=df, aes(x = ages, y = value, group = variable, 
                                 color = factor(variable)), size = lw, linetype = 2)
      }
    } # fit
    if (dta){
      ## add data
      if (ctr) {
        df = data.frame(ages = ages[,1], Ey = yy[,1,idx])
        df = melt(df,id="ages")
        # Create the plot
        if (dtatype!="l"){
        Plot <-Plot + 
          geom_point(data=df, aes(x = ages, y = value, group = variable, 
                                  color = factor(variable)), size = 1)
        } else {
          Plot <-Plot + 
            geom_line(data=df, aes(x = ages, y = value, group = variable, 
                                    color = factor(variable)), size = 0.5)
        }
      }
      
      # matlines(ages[,1], yy[,1,idx], col=1:I, lwd=0.53, lty=1,
      #                 type=dtatype,pch=20)
      if (case) {
        df = data.frame(ages = ages[,2], Ey = yy[,2,idx])
        df = melt(df,id="ages")
        # Create the plot
        if (dtatype!="l"){
          Plot <-Plot + 
            geom_point(data=df, aes(x = ages, y = value, group = variable, 
                                    color = factor(variable)), 
                       size = 2, shape=1)
        } else {
          Plot <-Plot + 
            geom_line(data=df, aes(x = ages, y = value, group = variable, 
                                    color = factor(variable)), 
                       size = 0.5, shape=1, linetype = 2) +
            theme(legend.position="none")
        }
      }
    } #dtad
  } # fit | dta
  if (prot){
    ## proteins w/o pat effects
    df= data.frame(ages, mp)
    Plot = ggplot(data=df, aes(x=X1,y=X1.1))+
      coord_cartesian(ylim=range(Eyp,na.rm=T), xlim=range(ages,na.rm=T))+
      ylab("mp") +xlab("ages")
    if (ctr) {
      df = data.frame(ages = ages[,1], Ey = Eyp0[,1,idx])
      df = melt(df,id="ages")
      # Create the plot
      Plot = Plot + 
        geom_line(data=df, aes(x = ages, y = value, group = variable, 
                               color = "grey"), size = 1) +
        theme(legend.position="none")+
        coord_cartesian(ylim=range(df$value,na.rm=T), xlim=range(ages,na.rm=T))
    }
    
    if (case) {
      df = data.frame(ages = ages[,2], Ey = Eyp0[,2,idx])
      df = melt(df,id="ages")
      # Create the plot
      Plot <-Plot + 
        geom_line(data=df, aes(x = ages, y = value, group = variable, 
                               color = "pink"), size = 1, linetype = 2)
    }
  }
  return(Plot)
}
