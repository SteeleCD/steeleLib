# run heirarchical clustering
runHclust = function(data,distMeth="euclidean",clustMeth="ward.D",doGroups=FALSE,nGroups=5,groups=NULL,file=NULL)
	{
	# distance and hclust
	distance = dist(data,method=distMeth)
	clusters = hclust(distance,method=clustMeth)
	# get groups
	if(doGroups)
		{
		groups = cutree(clusters,nGroups)
		}
	# get colour for label
	getCol = function(label,colGroups=groups,colours=c("green","red","black","blue","cyan"))
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
	return(clusters)
	}


