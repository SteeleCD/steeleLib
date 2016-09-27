# run heirarchical clustering
runHclust = function(data,distMethFun=NULL,distMeth="euclidean",clustMeth="ward.D2",doGroups=FALSE,nGroups=5,groups=NULL,file=NULL,plotcolours=c("green","red","black","blue","cyan"))
	{
	# distance and hclust
	if(is.null(distMethFun))
		{
		distance = dist(data,method=distMeth)
		} else {
		distance = distMethFun(data)
		}
	clusters = hclust(distance,method=clustMeth)
	# get groups
	if(doGroups)
		{
		groups = cutree(clusters,nGroups)
		groups = sapply(unique(groups),FUN=function(x) names(groups)[which(groups==x)])
		}
	# get colour for label
	getCol = function(label,colGroups=groups,colours=plotcolours)
		{
  		index = which(sapply(colGroups,FUN=function(x) label%in%x))
  		return(colours[index])
		}
	# colour the labels
	# adapted from:
	# https://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
	labelCol <- function(x) {
  		if (is.leaf(x)) {
    			## fetch label
    			label <- attr(x, "label") 
    			## set label color to red for A and B, to blue otherwise
    			attr(x, "nodePar") <- list(lab.col=sapply(label,getCol))
  			}
  		return(x)
		}
	# colour the labels
	dend = dendrapply(as.dendrogram(clusters), labelCol)
	# plot hclust
	if(!is.null(file)) pdf(file)
	plot(dend)
	if(!is.null(file)) dev.off()
	return(list(clusters=clusters,dend=dend))
	}

# K means AIC
kmeansAIC = function(fit)
	{
	m = ncol(fit$centers)
	n = length(fit$cluster)
	k = nrow(fit$centers)
	D = fit$tot.withinss
	return(data.frame(AIC = D + 2*m*k,
                  BIC = D + log(n)*m*k))
	}

# hclust followed by k means
hclustKmeans = function(distance,method="ward.D2",nGroups=4)
	{
	# heirarchical clustering
	clusters = hclust(distance,method)
	# hclust groups
	groups = cutree(clusters,nGroups)
	# centers of heirarchical clusters
	clust.centers = aggregate(pca$x,list(groups),mean)[,-1]
	# k means starting from HC centers
	KM = kmeans(pca$x,centers=clust.centers)
	return(KM)
	}



