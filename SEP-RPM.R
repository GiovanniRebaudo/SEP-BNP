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

# Run analysis, save results and produce some plots
# ex()

############ exploratory plots ###########
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
data_label2[data_country=='AFR'] = 'RU'
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

# Plot
w  <- read.myfile("w.txt",K,L) # w[k, iter, l]
niter  <- dim(w)[2]
pi <- read.myfile("pi.txt",K)
mu  <- read.myfile("mu.txt",L)
sig2  <- read.myfile("sig.txt",L)
Sj  <-  read.myfile("Sj.txt",K)
