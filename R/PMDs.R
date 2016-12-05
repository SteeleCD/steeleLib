# ====================================================================
#		Pipeline functions
# ====================================================================
# split betas by chromosome
splitBetas = function(chmrs,betas=betasSarc,intergenic=TRUE,mani)
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



