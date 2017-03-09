# run PC analysis
runPCA = function(file=NULL,data=NULL,fileOut=NULL)
	{
	if(is.null(file)&is.null(data)) stop("Please specify a data object or file")
	if(!is.null(file)&!is.null(data)) stop("Only specify a data object or a file, not both")
	if(!is.null(file)) 
		{
		load(file)
		data = norm$beta
		}
	pca = prcomp(t(data),retx=TRUE,center=TRUE,scale.=TRUE)
	pca$varianceExplained = (pca$sdev)^2 / sum(pca$sdev^2) 
	if(!is.null(fileOut))
		{
		save(list="pca",file=fileOut)
		}
	return(pca)
	}

# plot PCA results
plotPCA = function(pca,fileName=NULL,outDir,groups,
		colours=NULL,colFun=NULL,subset=NULL,
		arrows=TRUE,nArrows=10,xlim=NULL,
		ylim=NULL,PCs=c(1,2),amplifyArrows=4,
		legendloc="topleft")
	{
	if(!is.null(subset))
		{
		pca$x = pca$x[subset,]
		groups = groups[subset]
		}
	groups = as.factor(paste0(groups))
	if(is.null(colours))
		{
		if(is.null(colFun))
			{
			colours=rainbow(length(levels(groups)))
			} else {
			colours = do.call(colFun,list(n=length(levels(groups))))
			}
	}
	if(is.null(xlim)) xlim=range(pca$x[,PCs[1]])
	if(is.null(ylim)) ylim=range(pca$x[,PCs[2]])
	if(!is.null(fileName)) pdf(paste0(outDir,"/",fileName))
	plot(pca$x[,c(PCs[1],PCs[2])],col=colours[groups],xlim=xlim,ylim=ylim)
	if(arrows)
	{
	  arrowData = pcaArrows(pca,x=paste0("PC",PCs[1]),y=paste0("PC",PCs[2]),ntop=nArrows)
	  arrows(x0=rep(0,times=nArrows),y0=rep(0,times=nArrows),x1=arrowData$x*amplifyArrows,y1=arrowData$y*amplifyArrows)
	}
	legend(legendloc, legend = levels(groups), col = colours, cex = 0.8, pch = 1)
	if(!is.null(fileName)) dev.off()
}

# add arrows to PCA plot
pcaArrows = function(pca,x="PC1",y="PC2",ntop=10)
{
  data = data.frame(obsnames=row.names(pca$x),pca$x)
  datapc = data.frame(varnames=rownames(pca$rotation),pca$rotation)
  mult = min(
    (max(data[,y])-min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x])-min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc = transform(datapc,
                     v1=.7*mult*(get(x)),
                     v2=.7*mult*(get(y))
                     )
  lengths = sqrt(abs(datapc$v1)^2+abs(datapc$v2)^2)
  index = which(order(lengths)%in%c(1:ntop))
  print(cbind(paste0(datapc$varnames[index]),datapc$v1[index],y=datapc$v2[index]))
  return(list(var=paste0(datapc$varnames[index]),x=datapc$v1[index],y=datapc$v2[index]))
}

