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

# simulate new data from cox model and plot
plotCox = function(coxPH,data,whichPlot,whichKeep=NULL,outFile=NULL,colour=NULL)
	{
	# how many levels
	facLevels = levels(data[,whichPlot])
	nRep = length(facLevels)
	# get new data frame, keeping other variables constant
	#newData = sapply(data[,-grep(whichPlot,colnames(data))],FUN=function(x)
	#	if(is.factor(x)) return(rep(levels(x)[1],nRep)) else return(rep(mean(x,na.rm=TRUE),nRep))
	#	)
	#newData = cbind(newData,facLevels)
	#colnames(newData)[ncol(newData)] = whichPlot
	#newData = as.data.frame(newData)
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


