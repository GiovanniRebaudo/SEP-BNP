############ exploratory plots ###########

########## run after the microbiome-example.R
########## this is just for viewing the data

###### plot cdfs/pdfs of microbiomes in different country subjects
y=data3[,2:39]
gammas=colMeans(y)
y=sapply(1:n,function(i) as.matrix(y)[,i]/gammas[i])

y1=y[order(rowSums(y),decreasing = T),]

data_country=pivot_wider(data[,c(1:2,6)],names_from = Sample,values_from = nationality)
data_country=data_country[1,2:dim(data_country)[2]]

data_country=as.matrix(data_country)[1,]
names(data_country)=NULL

data_label=data_country
data_label[data_country=='AFR']=1
data_label[data_country=='AAM']=2

data_label2=data_country
data_label2[data_country=='AFR']='R-AF'
data_label2[data_country=='AAM']='AM-AF'

y_af=as.matrix(y1[,data_country=='AFR'])
y_am=as.matrix(y1[,data_country=='AAM'])


pdf_af=rowSums(y_af)/sum(rowSums(y_af))
pdf_am=rowSums(y_am)/sum(rowSums(y_am))

pdf_all=rowSums(y1)/sum(rowSums(y1))


ggplot()+
  geom_line(aes(0:119,c(0,cumsum(pdf_af)),col='R-AF'))+
  geom_line(aes(0:119,c(0,cumsum(pdf_am)),col='AM-AF'))+
  geom_line(aes(0:119,c(0,cumsum(pdf_all)),col='ALL'))+
  labs(x='OTU',y="Cumulative Relative Frequencies")+
  theme(legend.position = 'bottom')
ggsave("empiricalclustersgg.pdf",width=5,height=3,units = 'in')


dfplot_af=cbind(1:119,rowMeans(y_af)) %>% data.frame()
colnames(dfplot_af) =  c("OTU","Abundace") 

ggplot(dfplot_af,aes(OTU,Abundace))+
  geom_bar(stat='identity')+
  labs(x='OTU',y='Scaled Abundance')
ggsave("normedcount-hist-af.pdf",width=5,height=3.5,units = 'in')

dfplot_am=cbind(1:119,rowMeans(y_am)) %>% data.frame()
colnames(dfplot_am) =  c("OTU","Abundace") 

ggplot(dfplot_am,aes(OTU,Abundace))+
  geom_bar(stat='identity')+
  labs(x='OTU',y='Scaled Abundance')
ggsave("normedcount-hist-am.pdf",width=5,height=3.5,units = 'in')


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
tree_labels<- dendro_data(dend_expr, type = "rectangle")
tree_labels$labels <- cbind(tree_labels$labels, Subject = as.factor(data_label2))

ggplot() +
  geom_segment(data = segment(tree_labels), aes(x=x, y=y, xend=xend, yend=yend))+
  geom_segment(data = tree_labels$segments %>%
                 filter(yend == 0) %>%
                 left_join(tree_labels$labels, by = "x"), aes(x=x, y=y.x, xend=xend, yend=yend, color = Subject)) +
  #geom_text(data = label(tree_labels), aes(x=x, y=y, label=label, colour = Subject, hjust=-1), size=2) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  scale_colour_brewer(palette = "Dark2") + 
  theme_dendro() #+
  #ggtitle("Mayo Cohort: Hierarchical Clustering of Patients Colored by Diagnosis")


ggplot() +
  geom_segment(data = segment(tree_labels), aes(x=x, y=y, xend=xend, yend=yend))+
  geom_segment(data = tree_labels$segments %>%
                 filter(yend == 0) %>%
                 left_join(tree_labels$labels, by = "x"), aes(x=x, y=y.x, xend=xend, yend=yend)) +
  geom_text(data = label(tree_labels), aes(x=x, y=y, label=label, hjust=-1), size=2) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  scale_colour_brewer(palette = "Dark2") + 
  theme_dendro() #+
#ggtitle("Mayo Cohort: Hierarchical Clustering of Patients Colored by Diagnosis")


library("ape")
# Default plot
plot(as.phylo(hc2), cex = 0.6, label.offset = 0.5)
legend("topright",c(1,2),legend("R-AF","Af-AM"))

library(ggdendro)
ggdendrogram(hc2,cex=0.5)
ggdendrogram(hc2, rotate = TRUE, theme_dendro = T)
