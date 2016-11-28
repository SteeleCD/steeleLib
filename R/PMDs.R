# identify PMDs
library(steeleLib)
dirs = assignFolders()
# load betas
load(paste0(dirs$grpDir,"methylation/results/minfi/undiffSarc/norm/undiffNormOnlyNoFusions-betas-funnormSex-dropSex-dropLowQ.Rdata"))
# load manifest
manifest = read.csv(paste0(dirs$grpDir,"methylation/manifests/MethylationEPIC_v-1-0_B2.csv"),as.is=TRUE)

# ====================================================================
#		Pipeline functions
# ====================================================================
# split betas by chromosome
splitBetas = function(chmrs,betas=betasSarc,intergenic=TRUE,mani=manifest)
	{
	# probes that match chms
	index = which(mani$CHR==chmrs)
	subManifest = mani[index,]
	# intergenic
	if(intergenic)
		{
		index = which(subManifest$UCSC_RefGene_Name=="")
		subManifest = subManifest[index,]
		}
	# probes in array
	index = which(subManifest$IlmnID%in%rownames(betas))
	subManifest = subManifest[index,]
	# subset array
	subBetas = betas[subManifest$IlmnID,]
	# order
	subBetas = subBetas[order(subManifest$MAPINFO),]
	rownames(subBetas) = subManifest$MAPINFO[order(subManifest$MAPINFO)]
	return(subBetas)
	}

# get info for a window
getMeanSDn = function(start,end,pos,betas)
	{
	# which probes are in window
	index = which(pos>=start&pos<end)
	betas = betas[index]
	# mean, SD, n
	return(list(mean=mean(betas),sd=sd(betas),n=length(betas)))
	}

# moving window mean and SD
windowMeanSD = function(pos,betas,windowSize=50000)
	{
	# make windows
	pos = as.numeric(pos)
	starts = seq(from=min(pos),to=max(pos),by=windowSize)
	ends = starts+windowSize
	# get window mean etc
	means = mapply(FUN=getMeanSDn,start=starts,end=ends,MoreArgs=list(pos=pos,betas=betas))
	colnames(means) = starts
	return(means)
	}

# get CPs 
getCPs = function(betas,changeFilter=c(-0.1),colour="blue",
	windowSize=30000,methodCP="PELT",pen="MBIC",penVal=0,
	stat="Normal",plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotLines=TRUE,plotPoints=FALSE,plotAll=TRUE,nLimit=1,
	method="changepoint",segLimit=0)
	{
	library(changepoint)
	# split into p and q
	if(split)
		{
		centreIndex = which.max(diff(as.numeric(names(betas))))
		betas = list(betas[1:centreIndex],betas[(centreIndex+1):length(betas)])
		} else {
		betas = list(betas)
		}
	outAll = NULL
	for(i in 1:length(betas))
		{
		subBetas = betas[[i]]
		# get windowed means
		info = windowMeanSD(pos=names(subBetas),betas=subBetas,windowSize=windowSize)
		# remove NAs
		naIndex = !is.na(info["mean",])
		nIndex = info["n",]>=nLimit
		keepIndex = which(naIndex&nIndex)
		# get changepoints
		if(method=="changepoint")
			{
			# changepoints library
			CP = cpt.meanvar(as.numeric(info["mean",keepIndex]),
				Q=50,method=methodCP,penalty=pen,
				test.stat=stat,pen.value=penVal)
			# filters
			#meanB = mean(subBetas)
			index = table(c(which(diff(CP@param.est$mean)<changeFilter)#,# reduction > 0.1 
				#which(CP@param.est$mean<meanB) # seg mean less than average
				))
			index = as.numeric(names(index)[which(index==1)])
			# start and end points of CP
			starts = CP@cpts[index]
			ends = CP@cpts[index+1]
			# plot changepoints
			if(plotCP) 
				{
				plot(CP)
				abline(v=c(starts,ends),col="red",lty=2)
				}
			} else {
			# segment with DNAcopy
			seg = segment(CNA(logit(as.numeric(info["mean",keepIndex])),rep(1,times=length(keepIndex)),1:length(keepIndex)))
			# filter
			index = which(seg$output$seg.mean<segLimit)
			# start and end points of CP
			starts = seg$output[index,"loc.start"]
			ends = seg$output[index,"loc.end"]
			# plot changepoints
			if(plotCP) 
				{
				plot(seg)
				abline(v=c(starts,ends),col="red",lty=2)
				}
			}
		# genome scale start and end points
		startsGenome = as.numeric(colnames(info)[keepIndex][starts])
		endsGenome = as.numeric(colnames(info)[keepIndex][ends])	
		# plot changepoints
		if(plotGen)
			{
			plot(NA,xlim=range(as.numeric(names(subBetas))),ylim=0:1)
			mapply(FUN=function(a,b){polygon(x=c(a,a,b,b),y=c(-2,2,2,-1),col="gray")},
				a=startsGenome,b=endsGenome)
			points(names(subBetas),subBetas)
			lines(names(subBetas),smooth(subBetas),col=colour)
			}
		# return
		out = cbind(startsGenome,endsGenome)
		colnames(out) = c("start","end")
		outAll = rbind(outAll,out)
		}
	# plot combined
	if(plotAll)
		{
		betas = do.call(c,betas)
		plot(NA,xlim=range(as.numeric(names(betas))),ylim=0:1,ylab="Beta",xlab="Position")
		mapply(FUN=function(a,b){polygon(x=c(a,a,b,b),y=c(-2,2,2,-1),col="gray")},
			a=outAll[,1],b=outAll[,2])
		if(plotPoints) points(names(betas),betas)
		if(plotLines) lines(names(betas),smooth(betas),col=colour)
		}
	return(outAll)
	}

# ================================================================
# 			example
# ================================================================
betasAll = sapply(1:22,FUN=function(x) splitBetas(x,betas=betasSarc,mani=manifest,intergenic=FALSE),simplify=FALSE)
betasInter = sapply(1:22,FUN=function(x) splitBetas(x,betas=betasSarc,mani=manifest,intergenic=TRUE),simplify=FALSE)
betas = betasInter
# run CPs
par(mfrow=c(2,1))
METHOD="PELT"
WINDOW=30000
# 30000, 20000, 10000
change = c(-0.1)
change = 0
penalty="MBIC"
stat="Normal"
nLim=1
i=1
# changepoints
# tumour
getCPs(betas=betas[[i]][,1],col="red",
	windowSize=WINDOW,methodCP=METHOD,
	changeFilter=change,stat=stat,
	plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotPoints=TRUE,plotLines=FALSE,pen=penalty,
	penVal=0,plotAll=TRUE,nLimit=nLim)
# normal
getCPs(betas=betas[[i]][,ncol(betasSarc)],col="blue",
	windowSize=WINDOW,methodCP=METHOD,
	changeFilter=change,stat=stat,
	plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotPoints=TRUE,plotLines=FALSE,pen=penalty,
	penVal=0,plotAll=TRUE,nLimit=nLim)

# segment
# tumour
getCPs(betas=betas[[i]][,1],col="red",
	windowSize=WINDOW,methodCP=METHOD,
	changeFilter=change,stat=stat,
	plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotPoints=TRUE,plotLines=FALSE,pen=penalty,
	penVal=0,plotAll=TRUE,nLimit=nLim,method="segment")
# normal
getCPs(betas=betas[[i]][,ncol(betasSarc)],col="blue",
	windowSize=WINDOW,methodCP=METHOD,
	changeFilter=change,stat=stat,
	plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotPoints=TRUE,plotLines=FALSE,pen=penalty,
	penVal=0,plotAll=TRUE,nLimit=nLim,method="segment")

betas=betas[[i]][,ncol(betasSarc)];col="blue"
	windowSize=WINDOW;methodCP=METHOD
	changeFilter=change;stat=stat
	plotGen=FALSE;plotCP=FALSE;split=TRUE
	plotPoints=TRUE;plotLines=FALSE;pen=penalty
	penVal=0;plotAll=TRUE


getCPs = function(betas=betas,changeFilter=0,colour="blue",
	windowSize=30000,methodCP="PELT",pen="MBIC",
	stat="Normal",plotGen=FALSE,plotCP=FALSE,split=TRUE,
	plotLines=FALSE)

getCPs(betas[[1]][,1],plotCP=TRUE,changeFilter=c(-0.1),plotAll=TRUE)#,penVal=1,methodCP="BinSeg")
getCPs(betas[[5]][,"S00049501"],plotCP=TRUE,changeFilter=c(-0.1))



		CP = cpt.meanvar(as.numeric(info["mean",naIndex]),
			Q=50,method="BinSeg",penalty=pen,
			test.stat=stat,pen.value=0.1)
plot(CP)

# ================================================================
#		run pipeline
# ================================================================
betas = sapply(1:22,FUN=function(x) splitBetas(x,betas=betasSarc,mani=manifest,intergenic=TRUE),simplify=FALSE)
# get PMDs for each individual
method = "changepoint"
listAll = list()
for(i in 1:ncol(betasSarc))
	{
	listInd = list()
	#pdf(paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/all/",method,"/",colnames(betasSarc)[i],".pdf"))
	for(j in 1:22)
		{
		print(paste0(i,":",j))
		# compute actual PMDs
		#listInd[[j]] = getCPs(betas[[j]][,i],plotPoints=TRUE,plotLines=FALSE,method=method)
		listInd[[j]] = getCPs(betas[[j]][,i],plotAll=FALSE,method=method)
		}
	#dev.off()
	listAll[[i]] = listInd
	}
save(list="listAll",file=paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/PMDs-",method,".pdf"))
# count and plot PMDs
pca = prcomp(t(betasSarc),center=TRUE,scale=TRUE,retx=TRUE)
groups = cutree(hclust(dist(pca$x),method="ward.D2"),5)
counts = sapply(listAll,FUN=function(x) sapply(x, nrow))
pdf(paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/numDMPsPerGroup-",method,".pdf"))
boxplot(colSums(counts)~groups,border=1:5,ylab="# DMVs",xlab="Group")
dev.off()
pdf(paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/numDMPsPerChrom-",method,".pdf"))
barplot(rowSums(counts),names.arg=1:22,ylab="# DMVs",xlab="Chromosome")
dev.off()
# combine PMDs
all = NULL
for(i in 1:ncol(betasSarc))
	{
	for(j in 1:22)
		{
		tmp = listAll[[i]][[j]]
		if(nrow(tmp)!=0)
			{
			tmp = cbind(rep(colnames(betasSarc)[i],times=nrow(tmp)),
				rep(j,times=nrow(tmp)),
				tmp)
			all = rbind(all,tmp)
			}
		}
	}
colnames(all) = c("sample","chrom","start","end")
# get ready to abberation plot
library(copynumber)
forPlot = all
forPlot = cbind(all[,1:2],rep("p",nrow(all)),all[,3:4],rep(50,nrow(all)),rep(20,nrow(all)))
colnames(forPlot) = c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")
lims = sapply(betas,FUN=function(x) range(as.numeric(rownames(x))))
# fill missing segments
fillMissing = function(PMDs,start,end,sampleID,chrom)
	{
	if(nrow(PMDs)>0)
		{
		out = NULL
		for(i in 1:nrow(PMDs))
			{
			tmp = PMDs[i,,drop=FALSE]
			tmp[,"start.pos"] = tmp[,"end.pos"]
			if(i==nrow(PMDs))
				{
				tmp[,"end.pos"] = end
				} else {
				tmp[,"end.pos"] = PMDs[i+1,"start.pos"]
				}
			tmp[,"mean"] = 0
			out = rbind(out,PMDs[i,,drop=FALSE])
			out = rbind(out,tmp)		
			}
		tmp = PMDs[1,,drop=FALSE]
		tmp[,"end.pos"] = tmp[,"start.pos"]
		tmp[,"start.pos"] = start
		tmp[,"mean"] = 0
		out = rbind(tmp,out)
		} else {
		out = matrix(c(sampleID,chrom,"p",start,end,50,0),nrow=1)
		}
	colnames(out) = c("sampleID","chrom","arm","start.pos","end.pos","n.probes","mean")
	return(out)
	}
# fill the segments
allForPlot = sapply(colnames(betasSarc),FUN=function(x) sapply(1:22, FUN=function(y) {
PMDs = forPlot[which(forPlot[,"chrom"]==y&forPlot[,"sampleID"]==x),,drop=FALSE]
fillMissing(PMDs,lims[1,y],lims[2,y],x,y)},simplify=FALSE),simplify=FALSE)
allForPlot = do.call(rbind,sapply(allForPlot,FUN=function(x) do.call(rbind,x)))
# change to data frame
allForPlotDF = data.frame(sampleID=allForPlot[,"sampleID"],
	chrom=as.numeric(allForPlot[,"chrom"]),
	arm=allForPlot[,"arm"],
	start.pos=as.numeric(allForPlot[,"start.pos"]),
	end.pos=as.numeric(allForPlot[,"end.pos"]),
	n.probes=as.numeric(allForPlot[,"n.probes"]),
	mean=as.numeric(allForPlot[,"mean"]))
# finally, abberation plot
pdf(paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/aberrationPlot-",method,".pdf"))
plotAberration(allForPlotDF,thres.gain=0.3)
dev.off()
# aberration plot by group
for(i in 1:5)
	{
	pdf(paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/aberrationPlot-grp",i,"-",method,".pdf"))
	plotAberration(allForPlotDF[which(allForPlotDF$sampleID%in%names(groups)[which(groups==i)]),],thres.gain=0.3)
	dev.off()	
	}
# write GISTIC input
forGISTIC = allForPlotDF
forGISTIC = forGISTIC[,c("sampleID","chrom","start.pos","end.pos","n.probes","mean")]
colnames(forGISTIC) = c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
write.table(forGISTIC,paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/seg-",method,".txt"),sep="\t",row.names=FALSE,quote=FALSE)
# seg by group
for(i in 1:5)
	{
	write.table(forGISTIC[which(paste0(forGISTIC$Sample)%in%names(groups)[which(groups==i)]),],paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/seg",i,"-",method,".txt"),sep="\t",row.names=FALSE,quote=FALSE)	
	}
# =================================================================
# 			Development
# =================================================================
betas = sapply(1:22,splitBetas,simplify=FALSE)


pdf(paste0(dirs$dropDir,"PostDoc/methylation/plots/undiffSarc/bumps/PMDs/plotByChrom.pdf"))
par(mfrow=c(2,1),mar=c(5,4,0,0))
for(i in 1:22)
	{
	plot(rownames(betas[[i]]),betas[[i]][,75])
	#lines(rownames(betas[[i]]),smooth(betas[[i]][,75]),col="blue")
	plot(rownames(betas[[i]]),betas[[i]][,1])
	#lines(rownames(betas[[i]]),smooth(betas[[i]][,1]),col="red")
	}
dev.off()

library(MASS)
par(mfrow=c(2,1),mar=c(5,4,0,0))
dens = kde2d(x=as.numeric(rownames(betas[[i]])),y=betas[[i]][,75],n=300)
image(dens,col=c("purple","blue","white","orange","red"))
points(rownames(betas[[i]]),betas[[i]][,75])
dens = kde2d(x=as.numeric(rownames(betas[[i]])),y=betas[[i]][,1],n=300)
image(dens,col=c("purple","blue","white","orange","red"))
points(rownames(betas[[i]]),betas[[i]][,75])


infoTum = windowMeanSD(pos=rownames(betas[[1]]),betas=betas[[1]][,1])
par(mfrow=c(3,1),mar=c(5,4,0,0))
plot(rownames(betas[[i]]),smooth(betas[[i]][,75]),col="blue",type="l")
lines(rownames(betas[[i]]),smooth(betas[[i]][,1]),col="red")
plot(colnames(infoTum),infoTum["mean",])
plot(colnames(infoTum),infoTum["sd",])


infoNorm = windowMeanSD(pos=rownames(betas[[1]]),betas=betas[[1]][,75])
plot(rownames(betas[[i]]),smooth(betas[[i]][,75]),col="blue",type="l")
lines(rownames(betas[[i]]),smooth(betas[[i]][,1]),col="red")
plot(colnames(infoNorm),infoNorm["mean",])
plot(colnames(infoNorm),infoNorm["sd",])

index = which(!is.na(infoTum["mean",]))
plot(cpt.meanvar(as.numeric(infoTum["mean",index]),Q=50,method="PELT"))
index = which(!is.na(infoNorm["mean",]))
plot(cpt.meanvar(as.numeric(infoNorm["mean",index]),Q=50,method="PELT"))

index = which(!is.na(infoTum["mean",]))
CPtum = cpt.meanvar(as.numeric(infoTum["mean",index]),Q=50,method="PELT")
meanB = mean(betas[[1]][,1])
CPtum@param.est$mean

which(diff(CPtum@param.est$mean)<c(-0.1))
which(CPtum@param.est$mean<meanB)
index = table(c(which(diff(CPtum@param.est$mean)<c(-0.1)),# reduction > 0.1 
	which(CPtum@param.est$mean<meanB) # seg mean less than average
	))
index = as.numeric(names(index)[which(index==2)])
starts = CPtum@cpts[index]
ends = CPtum@cpts[index+1]
plot(CPtum)
abline(v=c(starts,ends),lty=2,col="red",lwd=2)

startsGenome = as.numeric(colnames(infoTum)[which(!is.na(infoTum["mean",]))][starts])
endsGenome = as.numeric(colnames(infoTum)[which(!is.na(infoTum["mean",]))][ends])
plot(NA,xlim=range(as.numeric(rownames(betas[[i]]))),ylim=0:1)
mapply(FUN=function(a,b){polygon(x=c(a,a,b,b),y=c(-2,2,2,-1),col="gray")},a=startsGenome,b=endsGenome)
lines(rownames(betas[[i]]),smooth(betas[[i]][,75]),col="blue")
lines(rownames(betas[[i]]),smooth(betas[[i]][,1]),col="red")


CPnorm = cpt.meanvar(as.numeric(infoNorm["mean",index]),Q=50,method="PELT")
meanB = mean(betas[[1]][,75])
CPnorm@param.est$mean


saveBetas = betas

start = min(as.numeric(rownames(betas[[1]])))
end = start+1000000
thisB = betas[[1]][which(as.numeric(rownames(betas[[1]]))<end&as.numeric(rownames(betas[[1]]))>=start),75]

start = 4.6e7
end = start+10000000
thisB = betas[[1]][which(as.numeric(rownames(betas[[1]]))<end&as.numeric(rownames(betas[[1]]))>=start),75]
plot(cpt.meanvar(thisB))



getCP = function(start,end,pos,betas)
	{
	index = which(pos>=start&pos<end)
	betas = betas[index]
	if(length(index)<20) return(NA)
	changepoint = cpt.meanvar(betas)
	return(changepoint)
	}

windowCP = function(pos,betas,windowSize=1000000)
	{
	pos = as.numeric(pos)
	starts = seq(from=min(pos),to=max(pos),by=windowSize)
	ends = starts+windowSize
	CPs = mapply(FUN=getCP,start=starts,end=ends,MoreArgs=list(pos=pos,betas=betas))
	return(CPs)
	}

cpTum = windowCP(pos=rownames(betas[[1]]),betas=betas[[1]][,1])

par(mfrow=c(1,length(cpTum)),mar=c(0,0,0,0))
sapply(1:length(cpTum),FUN=function(i) plot(cpTum[[i]],yaxt="n",ylab=NA))


par(mfrow=c(1,3),mar=c(0,0,0,0))
sapply(1:3,FUN=function(i) plot(cpTum[[i]],yaxt="n",ylab=NA))

par(mfrow=c(2,1))
plot(cpt.meanvar(betas[[1]][,1],method="BinSeg"))
plot(cpt.meanvar(betas[[1]][,75],method="BinSeg"))


par(mfrow=c(2,1))
plot(cpt.meanvar(betas[[1]][,1],method="PELT",Q=50))
plot(cpt.meanvar(betas[[1]][,75],method="PELT",Q=50))




     set.seed(25)
     genomdat <- rnorm(500, sd=0.1) +
     rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
     plot(genomdat)
     chrom <- rep(1:2,c(290,210))
     maploc <- c(1:290,1:210)
     test1 <- segment(CNA(genomdat, chrom, maploc))

set.seed(51)
     genomdat <- rnorm(500, sd=0.2) +
     rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),c(137,87,17,49,29,52,87,42))
     plot(genomdat)
     chrom <- rep(1:2,c(290,210))
     maploc <- c(1:290,1:210)
     test2 <- segment(CNA(genomdat, chrom, maploc))




rownames(manifest) = manifest$IlmnID
keptManifest = manifest[rownames(betasSarc),]
mySeg = segment(CNA(logit(betasSarc[,1]),keptManifest$CHR,keptManifest$MAPINFO))


info = windowMeanSD(pos=names(betas[[1]][,1]),betas=betas[[1]][,1],windowSize=30000)
naIndex = which(!is.na(info[1,]))
mySeg = segment(CNA(logit(as.numeric(info[1,naIndex])),rep(1,times=length(naIndex)),1:length(naIndex)))

info = windowMeanSD(pos=names(betas[[1]][,75]),betas=betas[[1]][,75],windowSize=30000)
naIndex = which(!is.na(info[1,]))
mySegNorm = segment(CNA(logit(as.numeric(info[1,naIndex])),rep(1,times=length(naIndex)),1:length(naIndex)))

