## Peak calling ##

#p_value = 0.05
#folder = 'AH6408C_MACS_wiggle_norm'

peak_call<-function(folder, p_value)
{
# load data (make sure the file with all chromosomes is not in the folder)

filenames = list.files(folder, full=T)[2:17]
alldata = lapply(filenames, read.table, skip=2, sep='\t')
names(alldata) = filenames
#for (i in 1:length(alldata))
#{
#  even = seq(nrow(alldata[[i]])) %% 2
#  alldata[[i]] = cbind(alldata[[i]][even==1,], alldata[[i]][even==0,])
#}


"
# if you want p_value all over the dataset
# sig
chr_pos=vector()
for (i in 1:16)
{
  chr_pos[i] = nrow(alldata[[i]])
}

allmatrix = alldata[[1]]
for (i in 1:15)
{
  allmatrix = rbind(allmatrix,alldata[[i+1]])
}

sig_pos = which(pnorm(allmatrix[,2], mean=mean(allmatrix[,2]),sd=sd(allmatrix[,2]),lower.tail=F)<0.15)

sig_data = list()
w=1
for (i in 1:16)
{
  chr_sig_index = sig_pos[which(sig_pos > w & sig_pos < chr_pos[i])]
  sig_data[[i]] = allmatrix[chr_sig_index,]
  w = chr_pos[i]
}
"

# call peaks
peaks = list()
for (i in 1:16)
{
  print(i)
#  if (nrow(sig_data[[i]])>0)
#  {
  sig_chr = alldata[[i]][which(pnorm(alldata[[i]][,2], mean=mean(alldata[[i]][,2]),sd=sd(alldata[[i]][,2]),lower.tail=F)<p_value),]
  index = vector()
  index[1] = 1
  k = 1
  for (j in 2:nrow(sig_chr))
  {
    if (sig_chr[j,1]-sig_chr[j-1,1] > 1)
    {
      k = k+1
      index[k] = j
    }
  }  
  new_matrix = data.frame(matrix(NA, ncol=5))
  colnames(new_matrix) = c('start','end','length','summit_position','summit_value')
  n = 0
  for (m in 1:(length(index)-1))
  {
    range = sig_chr[index[m]:(index[m+1]-1),]
    if (length(index[m]:(index[m+1]-1))>=150)
    {
      n = n+1
      new_matrix[n,1]=sig_chr[index[m],1]
      new_matrix[n,2]=sig_chr[index[m+1]-1,1]
      new_matrix[n,3]=new_matrix[n,2]-new_matrix[n,1]
      new_matrix[n,4]=mean(range[which(range[,2]==max(range[,2])),1])
      new_matrix[n,5]=max(range[,2])
    }
  }
  peaks[[i]]=new_matrix
#  }
}

# output a peaks.bed file
for (i in 1:16)
{
  if (i <= 9) {
    chr<-paste("chr0",i,sep="")
  } else {
    chr<-paste("chr",i,sep="")
  }
  if (is.null(peaks[[i]])==FALSE)
  {
    write.table(cbind(rep(chr,nrow(peaks[[i]])),peaks[[i]][,1:2]), file=paste(folder,'_peaks.bed',sep=""), quote=F, row.names=F, col.names=F, append=T, sep='\t')
  }
}
  

# output a summits.bed file
for (i in 1:16)
{
  if (i <= 9) {
    chr<-paste("chr0",i,sep="")
  } else {
    chr<-paste("chr",i,sep="")
  }
  if (is.null(peaks[[i]])==FALSE)
  {
   end = ceiling(peaks[[i]][,4])+1
   out = cbind(rep(chr,nrow(peaks[[i]])),ceiling(peaks[[i]][,4]))
   out = cbind(out,end)
   out = cbind(out,paste('peak',1:nrow(peaks[[i]]),sep=''))
   out = cbind(out,peaks[[i]][,5])
   write.table(out, file=paste(folder,'_summits.bed',sep=""), quote=F, row.names=F, col.names=F, append=T, sep='\t')
  }
}
}