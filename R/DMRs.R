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
	return(list(betas=betas,pos=info))
	}

plotDMR = function(betas,dmrs,index,manifest,flank=10000,groupIndices,doInvLogit=TRUE)
	{
	# get values in region
	region = getRegion(betas,dmrs[index,"chr"],dmrs[index,"start"],dmrs[index,"end"],manifest,flank)
	# get positions
	pos1 = matrix(region$pos,ncol=length(groupIndices[[1]]),nrow=nrow(region$betas))
	pos2 = matrix(region$pos,ncol=length(groupIndices[[2]]),nrow=nrow(region$betas))
	# inverse logit for Ms
	if(doInvLogit) region$betas = invlogit(region$betas)
	# plot
	plot(pos1,region$betas[,groupIndices[[1]]],ylim=range(region$betas),xlab="Position",ylab=ifelse(doInvLogit,"Beta","M"),main=paste0("Chromosome ",dmrs[index,"chr"]))
	points(pos2,region$betas[,groupIndices[[2]]],col="red")
	abline(v=c(dmrs[index,"end"],dmrs[index,"start"]),lty=2)
	}
