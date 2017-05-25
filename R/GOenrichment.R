# ================================================================
# GO enrichment from methylation data
# ================================================================

GOenrichmentMethy = function(geneFile=NULL,outDir,array="EPIC")
	{
	library(goseq)
	# get manifest
	manifest = steeleLib:::getManifestOld(array)
	# get all genes on array
	allGenes = sapply(manifest$UCSC_RefGene_Name,FUN=function(x) unique(strsplit(x,split=";")[[1]]))
	allGenes = table(unlist(allGenes))
	# read in DM genes
	DMG = read.table(geneFile,col.names="gene")
	if(nrow(DMG)==0) return()
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
	write.table(enrich[which(sigIndex&enrich$ontology=="BP"),],file=paste0(outDir,"/correctedSigGO-BP.txt"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrich[which(sigIndex&enrich$ontology=="CC"),],file=paste0(outDir,"/correctedSigGO-CC.txt"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrich[which(sigIndex&enrich$ontology=="MF"),],file=paste0(outDir,"/correctedSigGO-MF.txt"),row.names=FALSE,quote=FALSE,sep="\t")
	# no correction
	enrichNormal = goseq(pwf,"hg19","geneSymbol",method="Hypergeometric")
	sigIndexNorm = p.adjust(enrichNormal$over_represented_pvalue,method="BH")<=0.05
	write.table(enrichNormal[which(sigIndexNorm&enrichNormal$ontology=="BP"),],file=paste0(outDir,"/uncorrectedSigGO-BP.txt"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrichNormal[which(sigIndexNorm&enrichNormal$ontology=="CC"),],file=paste0(outDir,"/uncorrectedSigGO-CC.txt"),row.names=FALSE,quote=FALSE,sep="\t")
	write.table(enrichNormal[which(sigIndexNorm&enrichNormal$ontology=="MF"),],file=paste0(outDir,"/uncorrectedSigGO-MF.txt"),row.names=FALSE,quote=FALSE,sep="\t")
}

getGeneBackground = function(array)
{
  # get manifest
  manifest = steeleLib:::getManifestOld(array)
  # get all genes on array
  allGenes = sapply(manifest$UCSC_RefGene_Name,FUN=function(x) unique(strsplit(x,split=";")[[1]]))
  allGenes = table(unlist(allGenes))
  return(allGenes)
}

enrichment = function(genes,background)
{
  # get data in right format
  data = as.numeric(names(background)%in%genes)
  names(data) = names(allGenes)
  # p null
  pwf = nullp(genes,"hg19","geneSymbol",bias.data=as.integer(background))
  # enrichment
  enrich = goseq(pwf,"hg19","geneSymbol",method="Hypergeometric")
  # significantly enriched
  sigIndex = p.adjust(enrich$over_represented_pvalue,method="BH")<=0.05
  return(enrich[which(sigIndex),])
}


singleBootstrap = function(nGenes,geneBackground)
{
  genes = sample(names(geneBackground),nGenes)
  enrichment(genes,geneBackground)
}

bootstrapEnrichment = function(nGenes,geneBackground=NULL,nBootstrap=1000,array="EPIC")
{
  # gene background
  if(is.null(geneBackground)) geneBackground = getGeneBackground(array)
  # perform multiple enrichment analyses
  res = replicate(nBootstrap, singleBootstrap(nGenes,geneBackground))
}

GOenrichPipeline = function(geneFile=NULL,outDir,array="EPIC",nBootstrap=1000)
{
  library(goseq)
  # get background
  geneBackground = getGeneBackground(array)
  # read in DM genes
  DMG = read.table(geneFile,col.names="gene")
  if(nrow(DMG)==0) return()
  # enrichment
  enriched = enrichment(DMG$gene)
  # bootsrap
  if(nBootstrap>0)
  {
    bootstrap = bootstrapEnrichment(nGenes=nrow(DMG),
                                    geneBackground=geneBackground,
                                    nBootstrap=nBootstrap,
                                    array=array)
    enriched$FDR = sapply(enriched$category,FUN=function(x) 
          sum(bootstrap$category==x)/nBootstrap)
    return(enriched)
  } else {
    return(enriched)
  }
}
