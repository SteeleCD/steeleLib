# ==========================================================================
# functions for determing number of clusters using kMeans or hclust
# =========================================================================

# kMeans
fitKmeans = function(data,nProbes,PCH,main="test")
	{
	# get SD
	SDs = apply(data,MARGIN=1,sd)
	# n most variable probes
	indexRow = rev(order(SDs))[1:nProbes]
	# kmeans
	fits = sapply(2:10, FUN=function(x) kmeans(logit(t(combat$beta[indexRow,indexCol])),x))
	# AICs
	AICs = do.call(rbind,apply(fits,MARGIN=2,kmeansAIC))
	# choose minimum BIC
	fit = fits[,which.min(AICs[,2])]
	# pca
	pca = prcomp(logit(t(data[indexRow,])))
	# plot pca
	plot(pca$x[,1:2],col=fit$cluster,pch=as.numeric(PCH),main=main)
	# return clusters
	return(fit)
	}

# hclust
fitHclust = function(data,nProbes,PCH,main="test",distFun=function(x) cor(x),HCLUST="complete")
	{
	# get SD
	#SDs = apply(data,MARGIN=1,sd)
	#get mean absolute deviation
	SDs = apply(data,MARGIN=1,FUN=function(x) mean(abs(x-mean(x))))
	# n most variable probes
	indexRow = rev(order(SDs))[1:nProbes]
	# distance matrix
	dist = as.dist(distFun(data[indexRow,]))
	# heirarchical clustering
	clust = hclust(dist,method=HCLUST)
	fitInfo = sapply(2:10,FUN=function(x) cluster.stats(dist,cutree(clust,x)))
	cut = cutree(clust,(2:10)[which.max(unlist(fitInfo["pearsongamma",]))])
	# pca
	pca = prcomp(t(data[indexRow,]))
	# plot pca
	plot(pca$x[,1:2],col=cut,pch=as.numeric(PCH),main=main)
	# return clusters
	return(clust)
	}
