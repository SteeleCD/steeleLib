# Generate new data for cox, mean for numeric, single value for factors
genNewData = function(x,nRep,whichKeep=NULL)
	{
	if(is.factor(x)) 
		{
		if(is.null(whichKeep))
			{
			whichKeep=which.max(table(x))
			} else {
			whichKeep = which(levels(x)==whichKeep)
			}
		return(factor(rep(levels(x)[whichKeep],nRep),levels=levels(x)))
		} else {
		return(rep(mean(x,na.rm=TRUE),nRep))
		}
	}

# Generate new data for survival prediction
genNewDataAll = function(data,whichPlot,whichKeep=NULL)
	{
	# how many levels
	facLevels = levels(data[,whichPlot])
	nRep = length(facLevels)
	newData = data.frame(factor(facLevels,levels=facLevels))
	for(i in colnames(data))
		{
		if(i==whichPlot)
			{
			#newData[,i] = facLevels
			next
			} else {
			newData[,i] = genNewData(data[,i],nRep,whichKeep=whichKeep[[i]])
			}
		}
	colnames(newData)[1] = whichPlot
	return(newData)
	}

# simulate new data from cox model and plot
plotCox = function(coxPH,data,whichPlot,whichKeep=NULL,outFile=NULL,colour=NULL)
	{
	newData = genNewDataAll(data,whichPlot,whichKeep)
	# fit model again with new data
	fit <- survfit(coxPH, newdata = newData)
	# plot simulated data
	library(survminer)
	if(!is.null(outFile)) pdf(outFile)
	simSurv = ggsurvplot(fit, conf.int = FALSE, legend.labs=facLevels,
           ggtheme = theme_minimal(),palette=colour,linetype=1:nrow(newData))
		print(simSurv$plot)
	if(!is.null(outFile)) dev.off()
	return(newData)
	}



# simulate new data from AFT flexsurv model, and plot
plotAFT = function(object,whichPlot,data,times,whichKeep=NULL,log=TRUE,colour,legloc="bottomleft",...)
	{
	# simulate new data, then get survival probs
	info = summary(object,newdata=genNewDataAll(data,whichPlot,whichKeep=whichKeep))
	# get covariate names
	infoNames = sapply(names(info),FUN=function(x) strsplit(x,split=", ")[[1]],simplify=FALSE)
	infoNames = sapply(sapply(infoNames,FUN=function(x) x[grep(whichPlot,x)]),FUN=function(y) strsplit(y,split="=")[[1]][2])
	if(log)
		{
		# plot time on log scale
		plot(NA,xlim=range(log(times)),ylim=0:1,...)
		sapply(1:length(info),FUN=function(x) lines(info[[x]][,"time"],info[[x]][,"est"],col=colour[x],...))
		} else {
		# plot time on ordinary scale
		plot(NA,xlim=range(times),ylim=0:1,...)
		sapply(1:length(info),FUN=function(x) lines(exp(info[[x]][,"time"]),info[[x]][,"est"],col=colour[x],...))
		}
	# add legend
	legend(x=legloc,legend=infoNames,col=colour,lty=1,lwd=2)
	}

