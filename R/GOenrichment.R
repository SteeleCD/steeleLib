GOenrichmentMethy = function(geneFile=NULL,outDir,array="EPIC")
	{
	library(goseq)
	# get manifest
	manifest = steeleLib:::getManifest(array)
	# get all genes on array
	allGenes = sapply(manifest$UCSC_RefGene_Name,FUN=function(x) unique(strsplit(x,split=";")[[1]]))
	allGenes = table(unlist(allGenes))
	# read in DM genes
	DMG = read.table(geneFile)
	# put in right format for goseq
	data = as.numeric(names(allGenes)%in%DMG[,1])
	names(data) = names(allGenes)
	# generate null distribution
	pdf(paste0(outDir,"/p-fit.pdf"))
	pwf = nullp(data,"hg19","geneSymbol",bias.data=as.integer(allGenes))
	dev.off()
	# enrichment analysis
	enrich = goseq(pwf,"hg19","geneSymbol")
	# signifantly enriched
	sigIndex = p.adjust(enrich$over_represented_pvalue,method="BH")<0.05
	write.table(enrich[which(sigIndex&enrich$ontology=="BP"),],file=paste0(outDir,"/correctedSigGO-BP.csv"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrich[which(sigIndex&enrich$ontology=="CC"),],file=paste0(outDir,"/correctedSigGO-CC.csv"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrich[which(sigIndex&enrich$ontology=="MF"),],file=paste0(outDir,"/correctedSigGO-MF.csv"),row.names=FALSE,quote=FALSE,sep="\t")
	# no correction
	enrichNormal = goseq(pwf,"hg19","geneSymbol",method="Hypergeometric")
	sigIndexNorm = p.adjust(enrichNormal$over_represented_pvalue,method="BH")<=0.05
	write.table(enrichNormal[which(sigIndexNorm&enrichNormal$ontology=="BP"),],file=paste0(outDir,"/uncorrectedSigGO-BP.csv"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrichNormal[which(sigIndexNorm&enrichNormal$ontology=="CC"),],file=paste0(outDir,"/uncorrectedSigGO-CC.csv"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrichNormal[which(sigIndexNorm&enrichNormal$ontology=="MF"),],file=paste0(outDir,"/uncorrectedSigGO-MF.csv"),row.names=FALSE,quote=FALSE,sep="\t")
	}
