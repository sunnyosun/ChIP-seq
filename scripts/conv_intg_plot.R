# Signals at Convergent genes
# Sunny

# MNase data 3h
folder = 'AH119C_MACS_wiggle_norm/'
filenames = list.files(folder, full=T)[2:17]
alldata = lapply(filenames, read.table, sep='\t')
red1=alldata

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

## optional, getting the gene names for each cluster
#midmatrix=as.data.frame(matrix(NA,ncol=5,nrow=1))
#for (i in seq(1,(nrow(mid)-1),2))
#{
#  newline=c(mid[i,1],as.numeric(mid[i,5]),as.numeric(mid[i+1,4]),mid[i,10],mid[i+1,10])
#  midmatrix=rbind(midmatrix,newline)
#}
#midmatrix=midmatrix[2:nrow(midmatrix),]
#midmatrix[,6]=as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])
midmatrix=midmatrix[which(midmatrix[,5]>0),]

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
     if ((length(which(as.numeric(red1[[k]][,1])>=start & as.numeric(red1[[k]][,1])<=end)))!=0)
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
   bp=75
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
par(mfrow=c(4,4))
plot(alldata[[1]][,1]-1000,alldata[[1]][,2]/71.4, xlim=c(-500,500),pch=16,col=color[1],ylab='',xlab='chr1',type='b',ylim=c(0.8,1.4))
for (i in 2:16)
{
  plot(alldata[[i]][,1]-1000,alldata[[i]][,2]/71.4,col=color[i], xlim=c(-500,500),pch=16,type='b',xlab=paste('chr',i,sep=''),ylab='',ylim=c(0.8,1.5))
}
legend(locator(1),col=color,pch=16, legend=1:16, horiz=T, cex=0.85)

a=0
for (i in 1:16)
{
  b=mean(red1[[i]][,2])
  a=a+b
}
a/16=71.4

ave_mid=as.data.frame(matrix(NA,ncol=2,nrow=27))
for (i in 1:27)
{
  ave_mid[i,1]=alldata_mid[[1]][i,1]
  value=0
  for (k in c(1:16))
  {
    value=value+as.numeric(alldata_mid[[k]][i,2])
  }
  ave_mid[i,2]=value/16
}

par(mfrow=c(1,1))
plot(ave_mid[,1]-1000,ave_mid[,2]/71.4, xlim=c(-500,500),pch=16,col=color[1],ylab='Average MNase Signal',xlab='Distance from 3 ends of convergent genes',type='b',ylim=c(0.8,1.4))
