check = function(toCheck)
	{
	isNull = is.null(toCheck)
	if(isNull) return(TRUE)
	isEmpty = length(toCheck)==0
	if(isEmpty) return(TRUE)
	isNA = is.na(toCheck)
	if(isNA) return(TRUE)
	return(FALSE)
	}

combineSegs = function(chromData,nMin=15)
	{
	toSkip = NULL
	index = which(chromData[,5]<nMin)
	if(length(index)==0) return(chromData)
	for(i in 1:length(index))
		{
		distanceFromStart = abs(chromData[index[i],3]-chromData[index[i]-1,4])
		print(distanceFromStart)
		DS <<- distanceFromStart
		distanceFromEnd = abs(chromData[index[i],4]-chromData[index[i]+1,3])
		print(distanceFromEnd)
		DE <<- distanceFromEnd
		closest = which.min(c(distanceFromStart,distanceFromEnd))
		if(check(distanceFromStart)) closest=2
		if(check(distanceFromEnd)) closest=1
		if(check(distanceFromStart)&check(distanceFromEnd)) 
			{
			toSkip = c(toSkip,i)
			next
			}
		print(closest)
		CLOSEST  <<- closest
		if(closest==1) 
			{
			nprobes1 = chromData[index[i]-1,5]
			nprobes2 = chromData[index[i],5]
			logR1 = chromData[index[i]-1,6]
			logR2 = chromData[index[i],6]
			} else if(closest==2) {
			nprobes1 = chromData[index[i]+1,5]
			nprobes2 = chromData[index[i],5]
			logR1 = chromData[index[i]+1,6]
			logR2 = chromData[index[i],6]
			}
			chromData[index[i]-1,6]=(nprobes1*logR1+nprobes2*logR2)/(nprobes1+nprobes2)
			chromData[index[i]-1,5]=nprobes1+nprobes2
		}
	if(!is.null(toSkip)) index = index[-toSkip]
	chromData[-index,]
	}

plotFun = function(seg.mean,num.mark,seg.start,seg.end,
		YLIM,colours)
	{
	#scatterplot3d(z=seg.mean,y=num.mark,x=seg.position)
	plot(NA,ylim=range(seg.mean),xlim=range(c(seg.start,seg.end)),col=colours[num.mark],xaxt="n")
	abline(h=0)
	sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours[num.mark][x],lwd=3))
	}

plotByChrom = function(segFile=NULL,segObj=NULL,dataDir,
		outDir,fileName='segsByChrom.pdf')
	{
	library(copynumber)
	if(!is.null(segFile)) data = read.table(paste0(dataDir,"/",segFile),head=TRUE,sep="\t")
	if(!is.null(segObj)) data = segObj
	# plot by chms
	chms = unique(data[,2])
	individual = unique(data[,1])
	# combine small segs
	data = do.call(rbind,sapply(individual,FUN=function(x) do.call(rbind,sapply(chms,FUN=function(y) combineSegs(data[which(data[,1]==x&data[,2]==y),]),simplify=FALSE)),simplify=FALSE))
	# colours
	colours = colorRampPalette(c("red","green"))(max(data[,5]))
	library(scatterplot3d)
	# plot
	pdf(paste0(outDir,'/',fileName))
	sapply(individual,FUN=function(x) {
		par(mfrow=c(4,6),mar=c(0,0,0,0))
		sapply(chms,FUN=function(y) {
			print(paste0(x,":",y))
			index = which(data[,1]==x&data[,2]==y)
			plotFun(data[index,6],data[index,5],data[index,3],data[index,4],range(data[,6]),colours)
		})})
	dev.off()
	data = cbind(data[,1:2],rep("q",times=nrow(data)),data[,3:ncol(data)])
	colnames(data) = c("sampleID","chrom","arm","start.pos","end.pos","n.probes","logR.mean")
	pdf(paste0(outDir,'/copynumber.pdf'))
	plotAberration(data,thres.gain=0.2)
	dev.off()
	}
	
plotAbber = function(segFile=NULL,segObj=NULL,dataDir,
		outDir=NULL,fileName="abberationPlot.pdf",
		HEAD=TRUE,logCN=TRUE,thresh=0.2,combine=TRUE,
		chrom=NULL)
	{
	library(copynumber)
	if(!is.null(segFile)) data = read.table(paste0(dataDir,"/",segFile),head=HEAD,sep="\t")
	if(!is.null(segObj)) data = segObj
	if(!logCN) data[,6] = log2(data[,6])-1
	# plot by chms
	chms = unique(data[,2])
	individual = unique(data[,1])
	# combine small segs
	if(combine)
		{
		data = do.call(rbind,sapply(individual,FUN=function(x) do.call(rbind,sapply(chms,FUN=function(y) combineSegs(data[which(data[,1]==x&data[,2]==y),]),simplify=FALSE)),simplify=FALSE))	
		}
	data = cbind(data[,1:2],rep("q",times=nrow(data)),data[,3:ncol(data)])
	colnames(data) = c("sampleID","chrom","arm","start.pos","end.pos","n.probes","logR.mean")
	if(!is.null(outDir)) pdf(paste0(outDir,'/',fileName))
	plotAberration(data,thres.gain=thresh,chrom=chrom)
	if(!is.null(outDir)) dev.off()
	}

