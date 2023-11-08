############ exploratory plots ###########

########## run after the protein-example.R 
########## this is just for viewing the data

library(ggplot2)

##### plot difference between patients and control along age 
yy=data.frame(y)
colnames(yy)=colnames(y)
yy1=yy[,sample(1:dim(yy)[2],20,replace = F)]
condition=c(rep(0,16),rep(1,16))

yy1=cbind(condition,age,yy1) 

diff=y[17:32,]-y[1:16,]
dp=diff[,sample(1:dim(yy)[2],20,replace = F)]
for(i in 1:20){
  while(sum(is.na(dp[,i]))>0){
    dp[,i]=diff[,sample(1:dim(yy)[2],1)]
  }
}

pdf("yp-yc-truth-2.pdf",width = 5,height=3)
par(mar=c(4,4,3,2))
plot(age[1:16],dp[,1],type='l',xlab='Age',ylab="Log Abundance",ylim=c(-8,8),cex.axis=1,cex.lab=1)
for(i in 2:19){
  lines(age[1:16],dp[,i],col=i)
}
dev.off()


##### plot randomly selected proteins 
idx_plot=sample(1:dim(yy)[2],10,replace = F)
for(i in 1:10){
  while(sum(is.na(yy[,idx_plot[i]]))>0){
    idx_plot[i]=sample(1:dim(yy)[2],1)
  }
}
#ggplot(data=yy1)+
#  geom_line(mapping = aes(x = age, y = P31942,color=condition))
yy1=cbind(condition,age,y[,idx_plot[1:10]]) 

pdf("summary-random10-protein3.pdf",width=6,height=3.5)
#par(xpd=T, mar=par()$mar+c(0,0,0,3))
par(mar=c(4,4,3,2))
plot(yy1[17:32,2],yy1[17:32,3],type='l',ylim=c(-6,10),col=1,xlab='Age',ylab='Log Abundace',bty='l')            
lines(yy1[1:16,2],yy1[1:16,3],lty=2,col=1)
for(j in 2:10){
  lines(yy1[17:32,2],yy1[17:32,j+2],type='l',col=j)  
  lines(yy1[1:16,2],yy1[1:16,j+2],lty=2,col=j)
}
grid(5, 5, lty=1,lwd = 0.5)
dev.off()


  
