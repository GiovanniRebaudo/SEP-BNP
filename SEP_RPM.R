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
library(ggpubr)
library(salso)
library(reshape2)


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
  
  n <- dim(Y0)[2]                ## n patients (J in the paper's notation)
  remove_sparse = F
  if(remove_sparse){
    drop = apply(Y0==0,1,sum)>n/4  ## drop all OTUS with more than n/4 0's
    Y  <- Y0[!drop,]               ## dropping 20 of the originally 119 OTU's
    ## now total # 0's is only 5.
  } else {
    Y  <- Y0
  }
  gammas=colMeans(Y)
  Y=sapply(1:n,function(i) as.matrix(Y)[,i]/gammas[i])
  y=log(Y+0.05)                  ## log transform (we only have 5 0's)
  return(list(y=y,Y=Y)) ## return both - the log transf "y" and raw "Y"
}
## create global variables for data & fied dimensions
out <- read.dta() ## maybe as.matrix(..)?? See if needed..
y <-   out$y    ## log transf
Y  <-  out$Y    ## abs scale

N <- dim(y)[1]*dim(y)[2]  # global var's (I*J in the paper)
n <- dim(y)[2] # # patients (J in the paper)
B <- dim(y)[1] # # OTU's (I in the paper)

##### max size of the subject (K) and OTU (L) partition
K = 10
L = 15
## hyperprior pars
alpha = 1      ## for GEM prior on pi
beta  = 1      ## GEM prior on w[k], k = 1, ..., K

## hyperpars (use "empirical Bayes" i.e., match sample moments)
mu0  = mean(as.matrix(y))                  ## mu_l ~ N(mu0, sig0)
sig0 = 3  ## SD!

a0 = 3                                    ## sig2_l ~ Ga(a0/2, b0/2)
b0 = var(c(y))/4

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
  verbose_step = max(round(1e4/20),1)
  startTime = Sys.time()
  ex(niter  = 1e4, #iteration MCMC
     niter0 = 1e3, # estimating Sj after niter0 iterations and stop updating Sj
     niter1 = 2e3  # just for plot 
  ) 
  timeRPM = difftime(Sys.time(), startTime, units=("secs"))[[1]]
}

# Load output
pi    <- read.myfile("Data-and-Results/pi.txt")
Sj    <- read.myfile("Data-and-Results/Sj.txt")
w     <- read.myfile("Data-and-Results/w.txt", K, L) # w[k, iter, l]
mki   <- read.myfile("Data-and-Results/mki.txt", K)
mu    <- read.myfile("Data-and-Results/mu.txt")
sig2  <- read.myfile("Data-and-Results/sig.txt")
# nk=sapply(1:K,function(k) sum(Sj==k))
nk    <- read.myfile("Data-and-Results/nk.txt")
# nl=sapply(1:L,function(l) sum(c(mki)==l))
nl    <- read.myfile("Data-and-Results/nl.txt")
ll    <- read.myfile("Data-and-Results/logl.txt")
nkp   <- sum(nk>0)
nlp   <- sum(nl>0)
# summ  <- c(iter, nkp, nlp, nk, nl, ll, pmki)
summ  <- read.myfile("Data-and-Results/iter.txt")

niter <- dim(w)[2]
point_SJ = salso::salso(Sj, nRuns = 100, maxZealousAttempts = 100, loss=VI())


#plot
# Distributional clusters by OTU
Dietswap_dataset <- readRDS("Data-and-Results/Dietswap_dataset.RDS")
data = Dietswap_dataset[Dietswap_dataset$timepoint==1,]
data2=pivot_wider(data[,1:3],names_from = Sample,values_from = Abundance)
na_in_rows=sapply(1:dim(data2)[1],function(i) sum(data2[i,]==0))
data3=data.frame(data2[na_in_rows!=38,]) #remove all na rows
data=data3[,2:39]


# This plot is made adapting the code in Fradenti Github rep CommonAtomModel
cluster_per_obs <- rep(point_SJ, rep(nrow(data), ncol(data)))
LL <- apply(data,2,function(x) cumsum(sort((x),decreasing = T)) /sum((x)) )
dim(LL)
LL2 <- reshape2::melt(LL)
LL2 <- cbind(LL2,cluster_per_obs)
LL2 <- as_tibble(LL2) %>% mutate(cluster_per_obs=paste0("",cluster_per_obs))
P = ggplot(LL2)+
  geom_hline(yintercept = c(0,1),lty=3)+
  geom_line(aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),
            lwd=.5, alpha=.2)+
  geom_point(aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),
             alpha=.9)+
  theme_bw()+scale_color_manual("Subject\nCluster",
                                values = c(2,4,1,6,5)) + ylim(0, 1)+
  xlab("Taxa sorted by count")+ylab("Cumulative Relative Frequncy")
# ggsave(plot=P, file="Image/Sjclustergg2.pdf", height = 2.5, width = 4)    

#### Plot heatmaps
point_SJ = salso::salso(Sj, nRuns = 100, maxZealousAttempts = 100, loss=VI())
if(length(unique(point_SJ))!=3){ 
  print("Adapt plot code for number of clusters different than 3")
  stop()}

# 
niter  = 1e4 #iteration MCMC
niter0 = 1e3 # estimating Sj after niter0 iterations and stop updating Sj
niter1 = 2e3 # just for plot 
mki_all <- read.myfile("Data-and-Results/mki.txt",K)
summ = read.myfile("Data-and-Results/iter.txt")
it = summ[,1]
col = brewer.pal(9, name="YlOrRd")
## col = viridis(10)

M  <- nrow(mki_all)
M1 = which(it >= niter1)[1]
mki  <- array(0,dim=c(K,B,length(M1:M)))
idx = 0
for(m in M1:M){
  idx = idx + 1
  for(b in 1:B){
    mki[,b,idx] = mki_all[m, (b-1)*K+(1:K)]
  }
}

PSM1 = psm(t(mki[1,,]))
PSM2 = psm(t(mki[2,,]))
PSM3 = psm(t(mki[3,,]))

# Reorder rows and columns (observations) of a dissimilarity matrix intra groups 
# and possibly reorder also the groups (batch of observations)
reorder_dismat_j <-  function(dismat, groups, order.groups=NULL){
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

Plot_heat <- function(dissimlar_stable = dissimlar_stable,
                      I          = B){
  dismat      = round(dissimlar_stable, 2)
  dismat      = reorder_dismat_j(dismat,groups=rep(1,I))
  plot_dismat = reshape2::melt(dismat)
  ggplot(data=plot_dismat, aes(x=factor(Var1), y=factor(Var2), fill=value)) + 
    geom_tile() + theme_bw()+ 
    scale_y_discrete(limits=floor(seq(I, 1)),
                     breaks = floor(seq(1, I, length.out = 9)), 
                     labels = floor(seq(1, I, length.out = 9))) +
    scale_x_discrete(breaks = floor(seq(1, I, length.out = 9)), 
                     labels = floor(seq(1, I, length.out = 9))) +
    xlab("OTU") + ylab("OTU") +
    scale_fill_gradientn(colours = c("white", "yellow", "red"), 
                         values = rescale(c(0, 0.5, 1)), 
                         space = "Lab", name="") +
    theme(legend.position = "right", text = element_text(size=20))
}

P1 = Plot_heat(PSM1, B) + xlab("") + theme(legend.position="none")
P2 = Plot_heat(PSM2, B) + ylab("") + theme(legend.position="none")
P3 = Plot_heat(PSM3, B) + xlab("") + ylab("")



PHeat <- ggarrange(P1, P2, P3, labels = c("1", "2", "3"), nrow=1,
                   widths = c(1, 1, 1.2))

# ggsave(plot=PHeat, file="Image/mik_coclusterprob_sorted2.pdf", 
       # width=15, height=4.5)

Point_Mki1 = salso::salso(t(mki[1,,]), nRuns = 100, 
                          maxZealousAttempts = 100, loss=VI())
Point_Mki2 = salso::salso(t(mki[1,,]), nRuns = 100, 
                          maxZealousAttempts = 100, loss=VI())
Point_Mki3 = salso::salso(t(mki[1,,]), nRuns = 100, 
                          maxZealousAttempts = 100, loss=VI())

Zmelts <- rbind(
  reshape2::melt(y[order(Point_Mki1),point_SJ==1]) %>%
    mutate(point_SJ = "C-1"),
  reshape2::melt(y[order(Point_Mki2),point_SJ==2]) %>%
    mutate(point_SJ = "C-2") #,
  # cbind(Var1=1:length(Point_Mki3), Var2 = rep(1, length(Point_Mki3)), 
  #       reshape2::melt(Z[order(Point_Mki3),point_SJ==3])) %>%
  #   mutate(point_SJ = "C-3")
)
Zmelts$Var1 = as.factor(as.integer(Zmelts$Var1))

P = ggplot(Zmelts,aes(Var2,Var1,fill=value)) +
  geom_raster() + theme_bw()+ 
  facet_wrap(.~point_SJ, scales = "free")+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank())+
  labs(y= "OTU", x = "Subject") 

# ggsave(plot=P, file="Image/mb-heatmap-y2.pdf", height = 5, width = 6)    

