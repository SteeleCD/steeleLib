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
	return(list(betas=betas,pos=info,overlaps=overlapGenes))
	}

# function to plot DMRs with overlapping genes displayed (if any)
plotDMR = function(betas,dmrs,index,manifest,flank=10000,groupIndices,doInvLogit=TRUE)
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
	plot(pos1[toOrder,1],smooth(log(rowMeans(region$betas[,groupIndices[[1]]])[toOrder]/rowMeans(region$betas[,groupIndices[[2]]])[toOrder])),type="l",col="black",lwd=2,xlab=XLAB,ylab="log(ratio)",main=paste0("Chromosome ",dmrs[index,"chr"]),xaxt=XAXT)
	abline(h=0,lty=2)
	abline(v=c(dmrs[index,"end"],dmrs[index,"start"]),lty=2)
	if(length(region$overlaps$geneStart)>0)
		{
		par(mar=c(5,4,0,2))
		yVals = 0:length(region$overlaps$geneStart)
		plot(NA,xlim=range(pos1),ylim=range(yVals),yaxt="n",ylab=NA,xlab="Position")
		for(i in 1:length(yVals))
			{
			polygon(c(region$overlaps$geneStart[i],region$overlaps$geneStart[i],region$overlaps$geneEnd[i],region$overlaps$geneEnd[i]),c(yVals[i]+0.25,yVals[i]+0.75,yVals[i]+0.75,yVals[i]+0.25),col="black")
			text(x=region$overlaps$geneEnd[i]+500,y=yVals[i]+0.5,region$overlaps$geneName[i],cex=1.5)
			}
		}
	}
