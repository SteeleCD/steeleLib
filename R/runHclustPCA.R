# ==========================================================================
# Wrapper for FactoMineR
# ==========================================================================

library(FactoMineR)

runHCPC = function(pca,outDir,outFile,saveIntermediates=FALSE)
	{
	res = HCPC(as.data.frame(pca$x),nb.clust=-1,iter.max=10,min=3,max=NULL,graph=TRUE)
	if(saveIntermediates) save(list="res",file=paste0(outDir,"/HCPCres.Rdata"))
	pdf(paste0(outDir,"/",outFile))
	plot.HCPC(res,choice="tree")
	plot.HCPC(res,choice="bar")
	plot.HCPC(res,choice="map")
	plot.HCPC(res,choice="map",draw.tree=FALSE)
	plot.HCPC(res,choice="map",draw.tree=FALSE,ind.names=FALSE)
	plot.HCPC(res,choice="3D.map")
	dev.off()
	return(res)
	}


