setwd("~/MovingBlockBootstrap")
load("Data/MovingBlockBootstrap.RData")


#FORMATTING
dat = dat[order(dat$SiteNumber_new),]
Sites$names =Sites$SiteNumber_new
Sites=Sites[order(Sites$names),]
dat = cbind(dat,Sites[, c("x","y","names")])
library(mvabund)
set.seed(2477)
x=dat[,"x"]
y=dat[,"y"]
xy=data.frame(x,y)
dist.xy=dist(xy)
species=names(dat)[24:(length(names(dat))-3)]
dat1=dat[,species]
species=mvabund(dat1)


#MODEL: ADJUST FOR ENVIRONMENTAL VARIATION
model1=manyglm(species~poly(mint5,2)+poly(maxt95,2)+poly(AnnualMeanTemp,2)+poly(AnnualPrec,2),data=dat,family="binomial")
#model2=manyglm(species~poly(mint5,2)+poly(maxt95,2),data=dat,family="binomial")
#dim(residuals.manyglm(model1))
res=residuals.manyglm(model1)
#res2=residuals.manyglm(model2)

#Calculate PAIRWISE SITE CORRELATIONS
cor.site=matrix(nrow=nrow(dat1),ncol=nrow(dat1))
for( i in 1:3592){
  print(i)
  for (j in 1:i){
    cor.site[i,j]=
      cor(unlist(res[i,]),unlist(res[j,]))
    
  }
}

save(cor.site,file="cor.site.RData")
load("cor.site.RData")

#Order the pairwise distances
o=order(dist.xy)

#Order the pairwise correlations by the pairwise distances
cor.site.vec=c()
for (k in 1:3592){
  cor.site.vec=c(cor.site.vec,na.omit(cor.site[,k])[-1])
}
cor.site.ordered.by.dist=cor.site.vec[o]

dist.sorted=dist.xy[o]

#Create bins for the correlogram
#bins.1=c(0,200,400,800,1000,2000,4000,8000,16000,32000,64000,128000,301000)
bins.2=round(1000*seq(1,27,1)^1.75)

dist.bins2=cut(dist.sorted,bins.2,order=T)
levels(dist.bins2)=bins.2

dist.bins2.unfactor = as.numeric(as.character(dist.bins2))


#Calculate mean of pairwise correlations in each bin
cor.site.list2 = list()
n.in.bin=c()
mean_in_bin=c()
j=1
for (i in levels(dist.bins2)){
  which_in_bin = which(dist.bins2.unfactor == as.numeric(i))
  n.in.bin[j] = length(which_in_bin)
  mean_in_bin[j] = mean(cor.site.ordered.by.dist[which_in_bin])
  cor.site.list2 [[j]] = sum(n.in.bin*mean_in_bin)/sum(n.in.bin)
  j=j+1
}


#plot the average correlation of sites sepaarated by pairwise dist <d against d. 

plot(bins.2,unlist(cor.site.list2), ylab="mean correlation of site Dun Smyth resids", xlab="pairwise distance less than")

#on average, sites have average correlation of :
abline (h = mean(cor.site.ordered.by.dist), col="red")


add.random.line = function(){
cor.site.random = sample(cor.site.ordered.by.dist, size=length(cor.site.ordered.by.dist))
cor.site.list.random = list()
n.in.bin=c()
mean_in_bin=c()
j=1
for (i in levels(dist.bins2)){
  which_in_bin = which(dist.bins2.unfactor == as.numeric(i))
  n.in.bin[j] = length(which_in_bin)
  mean_in_bin[j] = mean(cor.site.random[which_in_bin])
  cor.site.list.random [[j]] = sum(n.in.bin*mean_in_bin)/sum(n.in.bin)
  j=j+1
}
lines(bins.2,unlist( cor.site.list.random), col="green")
;cor.site.list.random
}
set.seed(1052)
cor.site.list.random.store =list()
for(i in 1:100){
cor.site.list.random.store[[i]] =unlist(add.random.line())
}

y0 = min(sapply(cor.site.list.random.store,function(x){min(x[-length(x)])}))
y1= max(unlist(cor.site.list2)[-length(unlist(cor.site.list2))])
plot(bins.2,unlist(cor.site.list2), ylab="Average correlation of residuals", xlab="pairwise distance of sites (m)", ylim = c(y0,y1), type="o",cex.lab=1.5, cex.axis=1.5)
for(i in 1:100){
  lines(bins.2,unlist( cor.site.list.random.store[[i]]),col="grey")
}
abline(h=0, lwd=2)
abline (h = mean(cor.site.ordered.by.dist), col="red",lwd=2)
lines(bins.2,unlist(cor.site.list2), col="orange",type="o",lwd=2)


legend(100000,0.004,c("data", "average correlation" ,"correlation = 0","simulation envelope","(under independence)"), col = c("orange","red","black","grey",0), lwd =c(2,2,2,1,1),cex=1.5)
save.image(file="plotSiteCorrelation.RData")
