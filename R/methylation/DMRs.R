# function to get genes that overlap a region
getOverlapGenes = function(chr,start,end)
	{
	# get gRanges from this DMR
	ranges = gsub(" ","",paste0(chr,":",start,":",end))
	library(dplyr)
	gRanges = sapply(ranges, function (x) {res=strsplit(x, ':')}) %>%
		unlist %>%
		as.numeric %>%
		matrix(ncol=3, byrow=T) %>%
		as.data.frame %>%
		dplyr::select(chrom=V1, start=V2, end=V3) %>%
		mutate(chrom=paste0('chr', chrom)) %>%
		makeGRangesFromDataFrame
	# get Hsapiens genes
	library(Homo.sapiens)
	genesRanges = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	# overlaps between genes and DMR
	overlap = subsetByOverlaps(genesRanges,gRanges); overlap
	geneStart = overlap@ranges@start
	geneEnd = overlap@ranges@start+overlap@ranges@width
	OVERLAP <<- overlap
	# get gene names
	library(org.Hs.eg.db)
	INFO <<- org.Hs.egSYMBOL
	print(str(org.Hs.egSYMBOL))
	unmapped = org.Hs.eg.db::org.Hs.egSYMBOL
	mapped = mappedkeys(unmapped)
	genes = unlist(as.list(unmapped[mapped]))
	GENES <<- genes
	DMRgenes = genes[which(names(genes)%in%overlap$gene_id)]
	DMRGENES <<- DMRgenes
	# return
	return(list(geneStart=geneStart,geneEnd=geneEnd,geneName=DMRgenes))
	}

# function to get beta values of all probes within a region, and any genes that overlap the region
getRegion = function(betas,chr,start,end,manifest,flank=10000)
	{
	region = sort(c(start,end))
	region[1] = region[1]-flank 
	region[2] = region[2]+flank
	# down to chromosome
	manifest = manifest[which(manifest$CHR==chr),]
	# down to region
	index = which(manifest$MAPINFO>region[1]&manifest$MAPINFO<region[2])
	manifest = manifest[index,]
	# subset betas
	betas = betas[which(rownames(betas)%in%manifest$IlmnID),]
	# order
	betas = betas[order(rownames(betas)),]
	manifest = manifest[order(manifest$IlmnID),]
	# get info ready for output
	info = manifest$MAPINFO
	names(info) = manifest$IlmnID
	info = info[which(names(info)%in%rownames(betas))]
	# get overlapping genes
	overlapGenes = getOverlapGenes(chr,region[1],region[2])
	# return
	return(list(betas=betas,pos=info,overlaps=overlapGenes,manRegion=manifest))
	}

# function to plot DMRs with overlapping genes displayed (if any)
plotDMR = function(betas,dmrs,index,manifest,flank=10000,groupIndices,doInvLogit=TRUE,xOffset=10)
	{
	# get values in region
	region = steeleLib:::getRegion(betas,dmrs[index,"chr"],dmrs[index,"start"],dmrs[index,"end"],manifest,flank)
	# get positions
	pos1 = matrix(region$pos,ncol=length(groupIndices[[1]]),nrow=nrow(region$betas))
	pos2 = matrix(region$pos,ncol=length(groupIndices[[2]]),nrow=nrow(region$betas))
	# inverse logit for Ms
	if(doInvLogit) region$betas = invlogit(region$betas)
	# layout
	if(length(region$overlaps$geneStart)>0) 
		{
		layout(matrix(1:2,ncol=1,nrow=2),widths=1,heights=c(3,1))
		par(mar=c(0,4,4,2))
		XLAB = NA
		XAXT = 'n'
		} else {
		par(mar=c(5,4,4,2))
		layout(matrix(1,ncol=1,nrow=1),widths=1,heights=c(1))
		XLAB = "Position"
		XAXT = NULL
		}
	# plot DMR
	#plot(pos1,region$betas[,groupIndices[[1]]],ylim=range(region$betas),xlab=XLAB,ylab=ifelse(doInvLogit,"Beta","M"),main=paste0("Chromosome ",dmrs[index,"chr"]),col="blue",xaxt=XAXT)
	#points(pos2,region$betas[,groupIndices[[2]]],col="red")
	#abline(v=c(dmrs[index,"end"],dmrs[index,"start"]),lty=2)
	toOrder = order(pos1[,1])
	#lines(pos1[toOrder,1],rowMeans(region$betas[,groupIndices[[1]]])[toOrder],col="blue",lwd=2)
	#lines(pos2[toOrder,1],rowMeans(region$betas[,groupIndices[[2]]])[toOrder],col="red",lwd=2)
	# plot genes
	positions = pos1[toOrder,1]
	ratios = smooth(log(rowMeans(region$betas[,groupIndices[[1]]])[toOrder]/rowMeans(region$betas[,groupIndices[[2]]])[toOrder]))
	plot(positions,ratios,type="l",col="black",lwd=2,xlab=XLAB,ylab="log(ratio)",main=paste0("Chromosome ",dmrs[index,"chr"]),xaxt=XAXT)
	abline(h=0,lty=2)
	abline(v=c(dmrs[index,"end"],dmrs[index,"start"]),lty=2)
	if(length(region$overlaps$geneStart)>0)
		{
		par(mar=c(5,4,0,2))
		yVals = 1:length(region$overlaps$geneStart)
		YLIM = range(yVals)
		YLIM[1] = YLIM[1]-1
		yVals = yVals-1
		if(length(xOffset)==1) xOffset = rep(xOffset,times=length(region$overlaps$geneStart)) 
		plot(NA,xlim=range(pos1),ylim=YLIM,yaxt="n",ylab=NA,xlab="Position")
		for(i in 1:length(yVals))
			{
			xRange = range(positions)
			polygon(c(region$overlaps$geneStart[i],region$overlaps$geneStart[i],region$overlaps$geneEnd[i],region$overlaps$geneEnd[i]),c(yVals[i]+0.35,yVals[i]+0.65,yVals[i]+0.65,yVals[i]+0.35),col="black")
			if(region$overlaps$geneStart[i]<xRange[1]&region$overlaps$geneEnd[i]>xRange[2])
				{
				xText = mean(xRange)
				yText = yVals[i]+0.75
				} else if(region$overlaps$geneStart[i]<xRange[1]) {
				xText = region$overlaps$geneEnd[i]+abs(diff(xRange))/xOffset[i]
				yText = yVals[i]+0.5
				} else if(region$overlaps$geneEnd[i]>xRange[2]) {
				xText = region$overlaps$geneStart[i]-abs(diff(xRange))/xOffset[i]
				yText = yVals[i]+0.5
				} else {
				xText = region$overlaps$geneEnd[i]+abs(diff(xRange))/xOffset[i]
				yText = yVals[i]+0.5
				}
			text(x=xText,y=yText,region$overlaps$geneName[i])
			}
		}
	}

# volcano plot
plotVolcano = function(data,manifestProbes=NULL,file=NULL,PAR=list(mfrow=c(1,1),mar=c(5,4,4,2)),title="test",folder="/home/chris/Dropbox/PostDoc/methylation/plots/undiffSarc/DMP/volcano/",MAIN="test")
	{
	base = paste0(folder,title,"/")
	if(!is.null(file)) pdf(paste0(base,title,file))
	if(!is.null(PAR)) do.call(par,PAR)
	index = which(data$adj.P.Val<0.05)
	Xvals = data$logFC[index]
	Yvals = -log(data$adj.P.Val[index])
	probes = paste0(data$probeID[index])
	XLIM = range(Xvals)
	YLIM = range(Yvals)
	if(!is.null(manifestProbes))
		{
		manifestProbes = paste0(manifestProbes)
		Xvals = Xvals[which(probes%in%manifestProbes)]
		Yvals = Yvals[which(probes%in%manifestProbes)]
		probes = probes[which(probes%in%manifestProbes)]
		}
	plot(Xvals,Yvals,col=ifelse(abs(Xvals)<0.3,"black","red"),ylab="-log(P)",xlab="log(fold change)",xlim=XLIM,ylim=YLIM,main=MAIN)
	if(!is.null(file)) dev.off()
	return(list(hyper=probes[which(Xvals>0.3)],hypo=probes[which(Xvals<c(-0.3))]))
	}
