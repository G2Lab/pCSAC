#!/usr/bin/R
dens = 0.07
r0 = 60*dens
mu = r0+2.4
sigma = (r0+240)/6
source = "ENm008_GM12878"
source = "ENm008_K562"


SampleSize = 64
data <- read.table(paste(SampleSize,".comb.txt",sep=""))
nnode = max(data[,2])
nmax = max(data[,3])
freqMatrix = matrix(0,nrow=nnode,ncol=nnode)
for (i in 1:nrow(data)){
	v1 = data[i,1];
	v2 = data[i,2];
	v3 = data[i,3];
  freqMatrix[v1,v2] = freqMatrix[v2,v1] =  v3/nmax;
}
MeanRand = mean(freqMatrix)

# write.table(freqMatrix,paste(SampleSize,".freqMat.txt",sep=""),quote=F,row.names=F,col.names=F)
##########################################################################
# read length data
##########################################################################
lengthFile = "ENm008_GM12878_StartEnd.txt"
nodelength = read.table(lengthFile)
nodelength = (nodelength[,2]-nodelength[,1])*dens
#cat(nodelength)
plot(nodelength,main=paste(source,SampleSize ))
for(i in 1:length(nodelength)){
	lines(c(i,i),c(0,nodelength[i]),type="l")
}

##########################################################################
# read experimental data
##########################################################################
expdatafile = paste("../../ENm008/Nature_2010_ContactMap_",source,".txt",sep="")
expdata = read.table(expdatafile)
ProbExp = expdata/max(expdata)

nrow = nrow(expdata)
ncol = ncol(expdata)

PropensityMatrix = matrix(0,nrow=nrow,ncol=ncol)
for (i in 1:(nrow-2)){
	for (j in (i+2):nrow){
        if (freqMatrix[i,j]==0 && ProbExp[i,j] !=0){
            cat (paste(i,j, "freq=", freqMatrix[i,j], "Prob=", ProbExp[i,j],"\n"))
        }
		if (freqMatrix[i,j]==0){
			freqMatrix[i,j] = MeanRand
			#cat (paste(i,j,"\n"))
		}
		PropensityMatrix[i,j] = PropensityMatrix[j,i] = ProbExp[i,j]/freqMatrix[i,j]
	}
}


hist(PropensityMatrix[PropensityMatrix!=0],breaks=100,main=paste(source,SampleSize ))

PropensityMatrix = (PropensityMatrix/max(PropensityMatrix))*(2/sqrt(2*pi*sigma*sigma))

MaxP = max(PropensityMatrix)
MaxP

DistMatrix = matrix(0,nrow=nrow,ncol=ncol)
for (i in 1:(nrow-2)){
	for (j in (i+2):nrow){
		if (PropensityMatrix[i,j] != 0){
			# cat(paste(log( 
			# (2/sqrt(2*pi*sigma*sigma))*(MaxP/PropensityMatrix[i,j])),"\n"))

			DistMatrix[i,j] = DistMatrix[j,i] = mu + 
			sqrt(
				2*sigma*sigma*log( 
				(MaxP/PropensityMatrix[i,j])))
		}
	}
}

x = y = c(1:nrow(DistMatrix))
image(x,y,DistMatrix,main=paste(source,SampleSize ))
# 
hist(DistMatrix[DistMatrix!=0],breaks=100,main=paste(source,SampleSize))
text(50,70,paste("r0 =",r0))
text(50,65,paste("mu =",mu))
text(50,60,paste("sigma =",sigma))
# 
# edit(DistMatrix)
write.table(DistMatrix,paste(source,".dis.txt",sep=""),row.names = FALSE,
            col.names = FALSE)
#DistMatrix

# freqMatrix
