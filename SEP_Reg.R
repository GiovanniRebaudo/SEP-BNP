# Codes accompanying "Separate Exchangeability as 
# Modeling Principle in Bayesian Nonparametrics"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.15
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(tidyr)
library(mvtnorm)
library(invgamma)
library(splines)
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(reshape2)

# Data: 
#   is read in by readDta(), which is called right after it's defined. 
# Data comes from the file  "data_protein.RData"
# After calling readDta() the design matrix and responses are saved 
# in global variables "X, y"
# and C (=# columns=proteins), R (=# rows, i.e., patients) 
# and p (dim of design vector). 
# That is X is (R x p) etc.
#
#
# Prior: 
# prior par's are set in a global variable "prior = list(..)" 
# which is set right after reading in the data.
# In particular (m-eta, S-eta) are the normal moments of protein-specific 
# spline coefficients eta_j; (m-xi, S-xi) the same for patient-specific xi_i.
# (a-eta, b-eta) are the PY par's for proteins; 
# and (a-xi, b-xi) the same for patients.
# 
#
# main: 
# first sets up data structures for G_eta ~ PY(a,b) and G_xi ~ PY(a,b).
# Since we have y[ijt] ~ N(xi_i + x_i' eta_j, sig^2) 
# we introduce an adjusted outcome
# yt = y-xi with mdpOffset(mdpXi) before updating eta's, and similarly
# yt = y-eta with mdpOffset(mdpEta) before updating xi's.
# The yt is a global variable. 
# mdpInitClust(.) simply initializes the random partition of the eta_j's
# THe main MCMC is then
# for(it in 1:niter){
#   update eta
#   update xi
#   update sig^2
# }
#
#
# sim results: 
# are saved in mcmcUpd(). 
# All it does is just update running totals Ey (for fitted curves), 
# Eyp (for fitted curves *w/o* patient effects -- just proteins).
# Those are saved in "Ey.txt", "Ey2.txt"(for 2nd moments) etc.
# 
# 
# random partition: 
# mdp$s are random partition cluster membership indicators 
# (mdp=mdpEta, and mdp=mdpXi for proteins and patients, respectively).
# Those are saved in "sProt.txt" (just the first 250 proteins..), 
# and "sPat.txt", respectively.
# 
# 
# Polya urn: 
# see mdpUpdS(.) for the Polya urn for the PY processes.
# For clusters with cardinality < ns (=20) 
# I marginalize wrt cluster specific parameters. 
# That's a bit computation intensive (using candidate's formula). 
# For ns>20 i condition on cluster-specific pars. 
# Posterior uncertainty is  small it doesn't make a differcence for our data.
# It's a bit of a pain to keep track of marginal (or conditional, for ns>20) 
# cluster-specific likelihood under alternative cluster assignments. 
# To be checked further. Seems ok, also since the fitted curves seem ok.
# 
# 
# plots: 
# in the end you find some plotting routines. First call "pltInit(.)", 
# just to read in all the sim output from the files into global variables 
# (so we don't have to read in for each plot).
# plt_reg(.) is an omnibus plotting funciton. 
# Setting the arguments as desired you get whichever plots.
# See "paper_reg(.)" for what we did. Just ignore the "ps_reg(..)" funcitons. 
# That's my personal ones to get the right graphics par's
# maxDiff_reg() find the proteins with largest difference across time.

set.seed(1992)
source("SEP_fcts.R")

## global variables for the dta
out = readDta_reg(file="Data-and-Results/data_protein.RData") 
## make global vars for the data
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
    aeta = 1,        ## eta[j] ~ G_eta; G_eta ~ PY(aeta, N(meta,Seta))
    beta = 0.05, #- discount parameter
    ## row effects (patients)
    mxi = my, ## xi[i] ~ G_xi; G_xi ~ PY(axi, N(mxi,Sxi))
    Sxi = 25, # large to favor assignment of common effect to patients instead of protein
    axi = 0.1,
    bxi = -0.1, # # 0 as a DP
    ## hyperpars theta=sig2
    asig = 2,
    bsig = 2) ## w=1/sigs ~ Ga(asig/2, bsig/2); E(w) = m=a/b, V(w)=m/(b/2)



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

# Run MCMC
if (FALSE){
  startTime = Sys.time()
  main_reg(6000)
  timeREG = difftime(Sys.time(), startTime, units=("secs"))[[1]]
}

# Plots in the paper

pltInit_reg()
load("Data-and-Results/yt.RData")
of = maxDiff_reg()

P1 = plt_reg_ggplot(T,T,F,1,case=T,ctr=T)
P2 = plt_reg_ggplot(T,T,F,1001,case=T,ctr=T)

library(cowplot)

P1 = P1 + theme(axis.title = element_blank())
P2 = P2 + theme(axis.title = element_blank())
P  = plot_grid(P1, P2) 

# Individual plot       
P  = P +
  draw_label("ages", x= 0.52, y=  0, vjust=-0.5, angle= 0) +
  draw_label("Y",    x=  0, y=0.55, vjust= 1.5, angle=90)
P

#ggsave(plot=P, file ="Image/Ind_prot.pdf", 
#       width=20, height=8, units = 'cm')

######### Normal- QQtest #####################
chain = mcmc[,-1]
colnames(chain) = c("it", "SSM", "sig2", "K-prot", "K-pat", paste("nk",1:5,sep=""), paste("nk",1:5,sep=""))


# Ranking proteins

# Check
table(chain$`K-pat`)
table(chain$`K-prot`)


plt_reg(F,T,F,of[1:20],case=T,ctr=T,dtatype="l")
plt_reg_ggplot(F,T,F,of[1:20],case=T,ctr=T,dtatype="l")

summary(chain$`sig2`)


df <- data.frame(y = as.vector(yt)/mean(sqrt(chain$`sig2`)))
P <- ggplot( df, aes(sample = y))
P <- P + stat_qq() + stat_qq_line()

P

ggsave(plot=P, file ="Image/qq_prot.pdf", width=12, height=12, units = 'cm')

sProt <- read.csv("./Data-and-Results/sProt.txt", header=FALSE)
# Leggi senza virgole
colnames(sProt) = c("prot","cl")
sProt[,"prot"]  = factor(sProt[,"prot"]) 
Clust_prot_mat = matrix(NA, nrow=nrow(sProt)/C, ncol=C) 
colnames(Clust_prot_mat) = levels(sProt[,"prot"])
for (lev in levels(sProt[,"prot"])){
  Clust_prot_mat[,lev] = sProt[sProt[,"prot"]==lev,"cl"]
}


library(salso)
library(reshape2)  
library(scales)
library(plyr)          # version 1.8.8
# Compute probability of co-clustering of prot
dissimlar = psm(Clust_prot_mat)
VI_prot   = salso::salso(Clust_prot_mat,   loss=VI(), nCores = 4,  
                         maxNClusters = 21,
                         maxZealousAttempts = 100)
table(VI_prot)
dissimlar_ord_prot = reorder_dismat(dissimlar,VI_prot)


library(Cairo)
CairoPNG(filename ='Image/coclust_prot.png', width = 500, height = 500)
# Posterior probabilities of co-clustering of prot 
heatmap(dissimlar_ord_prot, Rowv = NA, Colv = NA, scale='none', 
        labRow = FALSE, labCol = FALSE,
        main="proteins prob coclust")
invisible(dev.off())

sPat <- read.csv("./Data-and-Results/sPat.txt", header=FALSE)
# Leggi senza virgole
colnames(sPat) = c("pat","cl")
sPat[,"pat"]  = factor(sPat[,"pat"]) 
Clust_pat_mat = matrix(NA, nrow=nrow(sPat)/R, ncol=R) 
colnames(Clust_pat_mat) = levels(sPat[,"pat"])
for (lev in levels(sPat[,"pat"])){
  Clust_pat_mat[,lev] = sPat[sPat[,"pat"]==lev,"cl"]
}

# Compute probability of co-clustering of pat
dissimlar = psm(Clust_pat_mat)
VI_pat    = salso::salso(Clust_pat_mat,   loss=VI(), nCores = 4,  
                         maxNClusters = 21,
                         maxZealousAttempts = 100)
table(VI_pat)
dissimlar_ord_pat  = reorder_dismat(dissimlar, VI_pat)

CairoPNG(filename ='Image/coclust_pat.png', width = 500, height = 500)
# Posterior probabilities of co-clustering of pat 
heatmap(dissimlar_ord_pat, Rowv = NA, Colv = NA, scale='none', 
        labRow = FALSE, labCol = FALSE,
        main="pat prob coclust")
invisible(dev.off())



