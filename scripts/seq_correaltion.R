# Correlation calculation between seq datasets OR chip datasets
# By Sunny

# loads the data for each chromosomes
# e.g. Smc4
folder = 'AH119C_MACS_wiggle_norm'
filenames = list.files(folder, full=T)[2:17]
alldata = lapply(filenames, read.table, skip=2, sep='\t')
names(alldata) = filenames
for (i in 1:length(alldata))
{
  even = seq(nrow(alldata[[i]])) %% 2
  alldata[[i]] = cbind(alldata[[i]][even==1,], alldata[[i]][even==0,])
}
AH6200C=alldata


# loads another datasets for each chromosomes
# e.g pH2A
folder = 'AH6408E_MACS_wiggle_norm'
filenames = list.files(folder, full=T)
alldata = lapply(filenames, read.table, skip=2, sep='\n')
names(alldata) = filenames
for (i in 1:length(alldata))
{
  even = seq(nrow(alldata[[i]])) %% 2
  alldata[[i]] = cbind(alldata[[i]][even==1,], alldata[[i]][even==0,])
}
# gives a name
ph2a=alldata


# column names 
# make sure the first columns are always called 'start', and the second columns are different between
# different samples
for (i in 1:16)
{
  colnames(AH119C[[i]])=c('start','119C')
}
for (i in 1:16)
{
  colnames(AH6200C[[i]])=c('start','6200C')
}


# if microarray data is from FE_K (lr)
dsb=list()
for (i in 1:16)
{
  dsb[[i]]=cbind(plat[chr(i),][chrorder(i),'start'],2^(lr[chr(i)][chrorder(i)]))
  colnames(dsb[[i]])=c('start','dsb')
}


# calculates the correlation coefficients between seq and chip data
for(i in 1:16)
{
  a<-merge(smc4[[i]],dsb[[i]])
  print (cor(a[,2],a[,3],use = "na.or.complete"))
}


# calculates the correlation coefficients between two seq datasets
c=vector()
for(i in 1:16)
{
  a<-merge(AH6200C[[i]],AH6200C_R[[i]])
  c[i]=(cor(a[,2],a[,3]))
}

# plotting
# plot
i=3
par(mfrow=c(3,1))
par(mar=c(0.5,4.1,1,2))
plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr[chr(i)][chrorder(i)]),frame.plot=F,pch=16, col='blue', yaxt='n',xaxt='n',ylab='',xlab='',xlim=c(min(plat[chr(i),][chrorder(i),'start']/1000),max(plat[chr(i),][chrorder(i),'start']/1000)),cex=1.5)
par(mar=c(0.5,4.1,1,2))
plot(red1[[i]][,1]/1000, red1[[i]][,2],frame.plot=F,pch=16, col='red', yaxt='n',xaxt='n',ylab='Signal',xlab='',xlim=c(min(plat[chr(i),][chrorder(i),'start']/1000),max(plat[chr(i),][chrorder(i),'start']/1000)),cex=0.8)
par(mar=c(5.5,4.1,1,2))
plot(smc4[[i]][,1]/1000, hop1[[i]][,2],frame.plot=F,pch=16, col='black', yaxt='n',ylab='',xlab=paste('Chromosome',i,' Position (kb)', sep=""),xlim=c(min(plat[chr(i),][chrorder(i),'start']/1000),max(plat[chr(i),][chrorder(i),'start']/1000)),cex=0.8)


# seq data
i=3
par(mfrow=c(2,1))
par(mar=c(1,4.1,5.5,2))
plot(rec114[[i]][,1]/1000, rec114[[i]][,2],frame.plot=F,pch=16, col='red',ylab='',cex=0.8, type='l',lwd=5, ylim=c(-0.5,5), xlab=paste('Chromosome',i,' Position (kb)', sep=""))
par(mar=c(5.5,4.1,1,2))
points(plat[chr(i),][chrorder(i),'start']/1000,2^(lr[chr(i)][chrorder(i)])-1,pch=16, col='black',ylab='',cex=0.8,lwd=5, ylim=c(0,5))
