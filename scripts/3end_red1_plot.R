# 3'-end regions analysis (Red1)
# Sunny

# load the gff file (saved as txt by excel)
out=read.table('SK1_annotation_modified.txt')

# extract the intergenic regions of convergent genes
newmatrix=matrix(NA,ncol=12,nrow=1)
for (i in 1:(nrow(out)-1))
{
  if (strsplit(as.character(out[i,10]),'')[[1]][1]=='Y' & strsplit(as.character(out[i+1,10]),'')[[1]][1]=='Y')
  {
    if (out[i,1]==out[i+1,1] & strsplit(as.character(out[i,10]),'')[[1]][7]=='W' & strsplit(as.character(out[i+1,10]),'')[[1]][7]=='C')
    {
      print (i)
      newline=out[i:(i+1),]
      newmatrix=rbind(newmatrix,newline)
    }
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
mid=newmatrix

# intergenic regions + 500bp on both sides
midmatrix=as.data.frame(matrix(NA,ncol=3,nrow=1))
for (i in seq(1,(nrow(mid)-1),2))
{
  newline=c(mid[i,1],as.numeric(mid[i,5])-500,as.numeric(mid[i+1,4])+500)
  midmatrix=rbind(midmatrix,newline)
}
midmatrix=midmatrix[2:nrow(midmatrix),]
midmatrix=cbind(midmatrix,2000/(as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])))
midmatrix[,5]=as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])-1000

# optional, getting the gene names for each cluster
midmatrix=as.data.frame(matrix(NA,ncol=5,nrow=1))
for (i in seq(1,(nrow(mid)-1),2))
{
  newline=c(mid[i,1],as.numeric(mid[i,5]),as.numeric(mid[i+1,4]),mid[i,10],mid[i+1,10])
  midmatrix=rbind(midmatrix,newline)
}
midmatrix=midmatrix[2:nrow(midmatrix),]
midmatrix[,6]=as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])
midmatrix=midmatrix[which(midmatrix[,6]>0),]

chr_len=read.table('sk1_chr.txt')
newmatrix=as.data.frame(matrix(NA,ncol=6,nrow=1))
colnames(newmatrix)=colnames(midmatrix)
for (i in 1:nrow(midmatrix))
{
  chr=midmatrix[i,1]
  len=chr_len[which(chr_len[,1]==chr),2]
  if (midmatrix[i,2]>10000 & midmatrix[i,3]<(len-10000))
  {
    newline=midmatrix[i,]
    newmatrix=rbind(newmatrix,newline)
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
midmatrix=newmatrix
midmatrix1=midmatrix

midmatrix_400_2000=midmatrix1[which(midmatrix1[,6]<2000 & midmatrix1[,6]>=400),]

midmatrix_250_400=midmatrix1[which(midmatrix1[,6]<400 & midmatrix1[,6]>=250),]

midmatrix_150_250=midmatrix1[which(midmatrix1[,6]<250 & midmatrix1[,6]>=150),]

midmatrix_100_150=midmatrix1[which(midmatrix1[,6]<150 & midmatrix1[,6]>=100),]

midmatrix_100=midmatrix1[which(midmatrix1[,6]<100),]

# remove those of negative lengths
midmatrix=midmatrix[which(midmatrix[,5]>0),]

chr_len=read.table('sk1_chr.txt')
newmatrix=as.data.frame(matrix(NA,ncol=5,nrow=1))
colnames(newmatrix)=colnames(midmatrix)
for (i in 1:nrow(midmatrix))
{
  chr=midmatrix[i,1]
  len=chr_len[which(chr_len[,1]==chr),2]
  if (midmatrix[i,2]>10000 & midmatrix[i,3]<(len-10000))
  {
    newline=midmatrix[i,]
    newmatrix=rbind(newmatrix,newline)
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
midmatrix=newmatrix
midmatrix1=midmatrix

# pick a threshold
midmatrix=midmatrix1[which(midmatrix1[,5]<2000 & midmatrix1[,5]>=400),]
t_400_2000=as.data.frame(table_400_2000)
t_400_2000[,2]=t_400_2000[,2]/sum(t_400_2000[,2])

midmatrix=midmatrix1[which(midmatrix1[,5]<400 & midmatrix1[,5]>=250),]
t_250_400=as.data.frame(table_250_400)
t_250_400[,2]=t_250_400[,2]/sum(t_250_400[,2])

midmatrix=midmatrix1[which(midmatrix1[,5]<250 & midmatrix1[,5]>=150),]
t_150_250=as.data.frame(table_150_250)
t_150_250[,2]=t_150_250[,2]/sum(t_150_250[,2])

midmatrix=midmatrix1[which(midmatrix1[,5]<150 & midmatrix1[,5]>=100),]
t_100_150=as.data.frame(table_100_150)
t_100_150[,2]=t_100_150[,2]/sum(t_100_150[,2])

midmatrix=midmatrix1[which(midmatrix1[,5]<100),]
t_100=as.data.frame(table_100)
t_100[,2]=t_100[,2]/sum(t_100[,2])

# plot average signals
alldata_mid=list()
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  chr=midmatrix[which(midmatrix[,1]==c),]
  if (nrow(chr)!=0)
    {wigmatrix=matrix(NA,ncol=2,nrow=1)
for (i in 1:nrow(chr))
{
  start=as.numeric(chr[i,2])
  end=as.numeric(chr[i,3])
  if ((length(which(red1[[k]][,1]>=start & red1[[k]][,1]<=end)))!=0)
  {
   tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,6]
   tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
   tmp=cbind(tmp1,tmp2)
   wigmatrix=rbind(wigmatrix,tmp)
  }
}
wigmatrix=wigmatrix[2:nrow(wigmatrix),]

# mean
wig_sort=wigmatrix[order(wigmatrix[,1]),]
wig_sort[,1]=round(wig_sort[,1])
mean_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2001))
for (i in 0:2000)
{
  mean_matrix[i+1,1]=i
  mean_matrix[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
}
#plot(mean_matrix, pch=16,cex=0.6,col='blue')

# smooth
bp=5
sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
for (i in seq(1,2000, bp))
{
  sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
  sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
}
alldata_mid[[k]]=sm_matrix
#plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
}
}
# plot each chromosome
alldata=alldata_mid
color_s=c('red','yellow','orange','blue','purple','springgreen','brown','burlywood','pink','grey','skyblue','turquoise','dark green','violetred','magenta','gold')
color=color_s
plot(alldata[[1]][,1]-1000,alldata[[1]][,2], pch=16,col=color[1],xlab='Distance from the midpoint of intergenic regions of convergent genes (bp)', ylab='Red1 Signal')
for (i in 2:16)
{
  points(alldata[[i]][,1]-1000,alldata[[i]][,2],col=color[i], pch=16)
}
legend(locator(1),col=color,pch=16, legend=1:16, horiz=T, cex=0.85)

# plot
i=1
par(mfrow=c(3,1))
par(mar=c(0.5,4.1,1,2))
plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr[chr(i)][chrorder(i)]),frame.plot=F,pch=16, col='blue', yaxt='n',xaxt='n',ylab='',xlab='',xlim=c(min(plat[chr(i),][chrorder(i),'start']/1000),max(plat[chr(i),][chrorder(i),'start']/1000)),cex=1.5)
par(mar=c(0.5,4.1,1,2))
plot(red1[[i]][,1]/1000, red1[[i]][,2],frame.plot=F,pch=16, col='red', yaxt='n',xaxt='n',ylab='Signal',xlab='',xlim=c(min(plat[chr(i),][chrorder(i),'start']/1000),max(plat[chr(i),][chrorder(i),'start']/1000)),cex=0.8)
par(mar=c(5.5,4.1,1,2))
plot(smc4[[i]][,1]/1000, smc4[[i]][,2],frame.plot=F,pch=16, col='purple', yaxt='n',ylab='',xlab=paste('Chromosome',i,' Position (kb)', sep=""),xlim=c(min(plat[chr(i),][chrorder(i),'start']/1000),max(plat[chr(i),][chrorder(i),'start']/1000)),cex=0.8)


# extract the intergenic regions of 5->3, 5->3 genes
newmatrix=matrix(NA,ncol=12,nrow=1)
for (i in 1:(nrow(out)-1))
{
  if (strsplit(as.character(out[i,10]),'')[[1]][1]=='Y' & strsplit(as.character(out[i+1,10]),'')[[1]][1]=='Y')
  {
    if (out[i,1]==out[i+1,1] & strsplit(as.character(out[i,10]),'')[[1]][7]=='W' & strsplit(as.character(out[i+1,10]),'')[[1]][7]=='W')
    {
      print (i)
      newline=out[i:(i+1),]
      newmatrix=rbind(newmatrix,newline)
    }
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
wa=newmatrix
wa=unique(wa)

# extract the intergenic regions of 3->5, 3->5
newmatrix=matrix(NA,ncol=12,nrow=1)
for (i in 1:(nrow(out)-1))
{
  if (strsplit(as.character(out[i,10]),'')[[1]][1]=='Y' & strsplit(as.character(out[i+1,10]),'')[[1]][1]=='Y')
  {
    if (out[i,1]==out[i+1,1] & strsplit(as.character(out[i,10]),'')[[1]][7]=='C' & strsplit(as.character(out[i+1,10]),'')[[1]][7]=='C')
    {
      print (i)
      newline=out[i:(i+1),]
      newmatrix=rbind(newmatrix,newline)
    }
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
cr=newmatrix
cr=unique(cr)


# intergenic regions + 500bp on both sides (watson)
wamatrix=as.data.frame(matrix(NA,ncol=3,nrow=1))
for (i in seq(1,(nrow(wa)-1),2))
{
  newline=c(wa[i,1],as.numeric(wa[i,5])-500,as.numeric(wa[i+1,4])+500)
  wamatrix=rbind(wamatrix,newline)
}
wamatrix=wamatrix[2:nrow(wamatrix),]
wamatrix=cbind(wamatrix,2000/(as.numeric(wamatrix[,3])-as.numeric(wamatrix[,2])))


# intergenic regions + 500bp on both sides (crick)
crmatrix=as.data.frame(matrix(NA,ncol=3,nrow=1))
for (i in seq(1,(nrow(cr)-1),2))
{
  newline=c(cr[i,1],as.numeric(cr[i,5])-500,as.numeric(cr[i+1,4])+500)
  crmatrix=rbind(crmatrix,newline)
}
crmatrix=crmatrix[2:nrow(crmatrix),]
crmatrix=cbind(crmatrix,2000/(as.numeric(crmatrix[,3])-as.numeric(crmatrix[,2])))


# plot wa
# plot average signals
alldata_wa=list()
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  chr=wamatrix[which(wamatrix[,1]==c),]
  wigmatrix=matrix(NA,ncol=2,nrow=1)
  for (i in 1:nrow(chr))
  {
    start=as.numeric(chr[i,2])
    end=as.numeric(chr[i,3])
    if ((length(which(red1[[k]][,1]>=start & red1[[k]][,1]<=end)))!=0)
    {
      tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,4]
      tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
      tmp=cbind(tmp1,tmp2)
      wigmatrix=rbind(wigmatrix,tmp)
    }
  }
  wigmatrix=wigmatrix[2:nrow(wigmatrix),]
  
  # mean
  wig_sort=wigmatrix[order(wigmatrix[,1]),]
  wig_sort[,1]=round(wig_sort[,1])
  mean_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2001))
  for (i in 0:2000)
  {
    mean_matrix[i+1,1]=i
    mean_matrix[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
  }
  #plot(mean_matrix, pch=16,cex=0.6,col='blue')
  
  # smooth
  bp=5
  sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
  for (i in seq(1,2000, bp))
  {
    sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
    sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
  }
  alldata_wa[[k]]=sm_matrix
  #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
}
alldata=alldata_wa
color_s=c('red','yellow','orange','blue','purple','springgreen','brown','burlywood','pink','grey','skyblue','turquoise','dark green','violetred','magenta','gold')
color=color_s
plot(alldata[[1]][,1]-1000,alldata[[1]][,2], pch=16,ylim=c(-1,6),col=color[1],xlab='Distance from the midpoint of intergenic regions of convergent genes (bp)', ylab='Red1 Signal')
for (i in 2:16)
{
  points(alldata[[i]][,1]-1000,alldata[[i]][,2],col=color[i], pch=16)
}
legend(locator(1),col=color,pch=16, legend=1:16, horiz=T, cex=0.85)


# plot cr
# plot average signals
alldata_cr=list()
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  chr=crmatrix[which(crmatrix[,1]==c),]
  wigmatrix=matrix(NA,ncol=2,nrow=1)
  for (i in 1:nrow(chr))
  {
    start=as.numeric(chr[i,2])
    end=as.numeric(chr[i,3])
    if ((length(which(red1[[k]][,1]>=start & red1[[k]][,1]<=end)))!=0)
    {
      tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,4]
      tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
      tmp=cbind(tmp1,tmp2)
      wigmatrix=rbind(wigmatrix,tmp)
    }
  }
  wigmatrix=wigmatrix[2:nrow(wigmatrix),]
  
  # mean
  wig_sort=wigmatrix[order(wigmatrix[,1]),]
  wig_sort[,1]=round(wig_sort[,1])
  mean_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2001))
  for (i in 0:2000)
  {
    mean_matrix[i+1,1]=i
    mean_matrix[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
  }
  #plot(mean_matrix, pch=16,cex=0.6,col='blue')
  
  # smooth
  bp=5
  sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
  for (i in seq(1,2000, bp))
  {
    sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
    sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
  }
  alldata_cr[[k]]=sm_matrix
  #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
}
alldata=alldata_cr
color_s=c('red','yellow','orange','blue','purple','springgreen','brown','burlywood','pink','grey','skyblue','turquoise','dark green','violetred','magenta','gold')
color=color_s
plot(alldata[[1]][,1]-1000,alldata[[1]][,2], pch=16,ylim=c(-1,6),col=color[1],xlab='Distance from the midpoint of intergenic regions of convergent genes (bp)', ylab='Red1 Signal')
for (i in 2:16)
{
  points(alldata[[i]][,1]-1000,alldata[[i]][,2],col=color[i], pch=16)
}
legend(locator(1),col=color,pch=16, legend=1:16, horiz=T, cex=0.85)


# average signal of all chromosomes
# convergent 3'-end genes
ave_mid=as.data.frame(matrix(NA,ncol=2,nrow=400))
for (i in 1:400)
{
  ave_mid[i,1]=alldata_mid[[1]][i,1]
  value=0
  for (k in c(1:16))
  {
    value=value+as.numeric(alldata_mid[[k]][i,2])
  }
  ave_mid[i,2]=value/16
}

ave_mid_150_250=ave_mid
ave_mid_100_150=ave_mid
ave_mid_100=ave_mid
ave_mid_250_400=ave_mid
ave_mid_400_2000=ave_mid

plot(ave_mid_150_250[,1],ave_mid_150_250[,2],col='red',frame.plot=F,pch=16,ylim=c(0,5),
     xlab='Intergenic region (bp)', ylab='Red1 Signal')
points(ave_mid_100_150[,1],ave_mid_100_150[,2],col='blue',pch=16)
points(ave_mid_100[,1],ave_mid_100[,2],col='orange',pch=16)
points(ave_mid_250_400[,1],ave_mid_250_400[,2],col='green',pch=16)
points(ave_mid_400_2000[,1],ave_mid_400_2000[,2],col='purple',pch=16)

legend(locator(1),pch=16,col=c('orange','blue','red','green','purple'),legend=c('<100: 137','100-150: 176','150-250: 242','250-400: 240','400-2000: 143'))


# watson
ave_wa=as.data.frame(matrix(NA,ncol=2,nrow=400))
for (i in 1:400)
{
  ave_wa[i,1]=alldata_wa[[1]][i,1]
  value=0
  for (k in 1:16)
  {
    value=value+alldata_wa[[k]][i,2]
  }
  ave_wa[i,2]=value/16
}

# crick
ave_cr=as.data.frame(matrix(NA,ncol=2,nrow=400))
for (i in 1:400)
{
  ave_cr[i,1]=alldata_cr[[1]][i,1]
  value=0
  for (k in 1:16)
  {
    value=value+alldata_cr[[k]][i,2]
  }
  ave_cr[i,2]=value/16
}

# watson+crick
ave_cr_r=cbind(ave_cr[,1],ave_cr[400:1,2])
wc=cbind(ave_wa[,1],(ave_wa[,2]+ave_cr_r[,2])/2)



# plot average signal of all chromosomes
# stop codon for mid is 611
# stop codon for wc is 444
# 167
mid_plot=ave_mid[which(ave_mid[,1]>=611),]
wc_plot=wc[which(wc[,1]>=444),]

plot(mid_plot[,1]-613,mid_plot[,2],col='red',frame.plot=F,pch=16,ylim=c(0,1.5),
     xlab='Distance of the stop codon away from the gene (bp)', ylab='Red1 Signal')
points(wc_plot[,1]-448,wc_plot[,2],col='blue',pch=16)
legend(locator(1),pch=16,col=c('red1','blue'),legend=c('convergent genes','unidirectional genes'))


midmatrix2=as.data.frame(matrix(NA,ncol=3,nrow=1))
colnames(midmatrix2)=colnames(midmatrix)
for (i in 1:nrow(midmatrix))
{
  if (midmatrix[i,2]<midmatrix[i,3])
  {
    midmatrix2=rbind(midmatrix2,midmatrix[i,])
  }
}
midmatrix2=midmatrix2[2:nrow(midmatrix2),]


