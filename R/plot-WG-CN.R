# =================================================================================
#			NECESSARY FOR NAMES
# =================================================================================
# 			CONVERT NAME TO HGNC
# ================================================================================
library(GenomeGraphs)
setMethod("initialize", "GeneRegion", function(.Object,...){
    .Object <- callNextMethod()
    strand  = switch(.Object@strand, "+" = 1, "-" = -1)
    
    ##-- changing start and end positions to capture genes on the edges. 
    .Object@start <- .Object@start - 2000
    .Object@end <- .Object@end + 2000

    if (!is.null(.Object@biomart)) {
       # .Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank", "structure_transcript_chrom_strand","structure_biotype"),filters=c("chromosome_name", "start", "end", "strand"), values=list(.Object@chromosome,.Object@start, .Object@end, strand), mart=.Object@biomart)
ens = getBM(c("external_gene_name","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank", "strand"),filters=c("chromosome_name", "start", "end", "strand"), values=list(.Object@chromosome,.Object@start, .Object@end, strand), mart=.Object@biomart)      
#tmp = ens[order(ens$exon_chrom_start),]
#TMP <<- tmp
#genes = unique(tmp$external_gene_name)
#GENES <<-genes
#for(i in genes) 
#	{
#	index = which(tmp$external_gene_name==i)
#	if(length(index)>1) 
#		{
#		tmp$external_gene_name[index[-round(median(1:length(index)))]] = ""
#		}
#	}
 .Object@ens <- ens#getBM(c("external_gene_name","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank", "strand"),filters=c("chromosome_name", "start", "end", "strand"), values=list(.Object@chromosome,.Object@start, .Object@end, strand), mart=.Object@biomart)
         if(!is.null(.Object@ens)){
     .Object@ens <- cbind(.Object@ens, biotype=rep("protein_coding", length(.Object@ens[,1])))
   }
    }
    
    if (is.null(.Object@ens)) {
        setPar(.Object, "size", 0)
    }
    .Object
})

makeGeneRegion <- function(start, end, chromosome, strand, biomart, dp = NULL){
    if(missing(start)) stop("Need to specify a start for creating a GeneRegion")
    pt <- getClass("GeneRegion")@prototype
    if (is.null(dp))
        dp <- pt@dp
    if(is.numeric(chromosome))
        chromosome = as.character(chromosome)
    new("GeneRegion", start = start, end = end, chromosome = chromosome, strand = strand ,biomart = biomart, dp = dp)
}

# ====================================================================================
#		# ONLY PRINT NAME ONCE FOR EACH GENE
# ====================================================================================
.drawGene <- function(gdObject, minBase, maxBase, vpPosition) {
    ens <- GenomeGraphs:::getExonModel(gdObject)
    if (dim(ens)[1] == 0) {
       # warning("No genes in gene region.")
        return(NULL)
    }
    pushViewport(dataViewport(xData=c(minBase, maxBase), extension = 0,
                              clip = TRUE, yscale = c(0, 40),
                              layout.pos.col=1, layout.pos.row = vpPosition))
    
    if(ens[1,6] == 1) { pp = 25 }
    else { pp = 15 }

    for(i in seq(along=ens[,1])) {
 
        color <- GenomeGraphs:::getColor(gdObject)
       
        if (!is.null(GenomeGraphs:::getBiotypeColor(gdObject, as.character(ens[i,8]))))
            color <- GenomeGraphs:::getBiotypeColor(gdObject, as.character(ens[i,8]))
        grid.rect(ens[i, 5], 5, width = ens[i, 5] - ens[i, 4], height = 30,
                  gp=gpar(col = "black", fill = color), default.units="native",
                  just = c("right", "bottom"))
     
      
    }

    
    genes = unique(ens[,1])
    for(g in seq(along=genes)){
        color = GenomeGraphs:::getColor(gdObject)
        exons = ens[ens[,1]==genes[g],-c(1,2,3,6)]
        ord = order(exons[,1])
        exons = exons[ord,]
        if(!is.null(GenomeGraphs:::getBiotypeColor(gdObject,as.character(exons[1,4])))) color = GenomeGraphs:::getBiotypeColor(gdObject,as.character(exons[1,4]))
        for(j in seq(along=exons[,1])){
            if(j < length(exons[,1])){
                grid.lines(c(exons[j,2],exons[j,2]+((exons[j+1,1] - exons[j,2])/4)),c(20,pp),
                           default.units = "native", gp=gpar(col = color, just = c("right", "bottom")))
                grid.lines(c(exons[j,2]+((exons[j+1,1] - exons[j,2])/4), exons[j+1,1]),c(pp,20),
                           default.units = "native", gp=gpar(col = color, just = c("right", "bottom")))
            }
        }
  if (GenomeGraphs:::getPlotId(gdObject)) {
            rot <- getPar(gdObject, "idRotation")
            col <- getPar(gdObject, "idColor")
GENES <<- genes
GENE <<- g
TMP <<- ens
index = which(ens[, 1]==genes[g])
starts = median(ens[index, 4])
ends = median(ens[index, 5])
            grid.text(genes[g], (starts + ends)/2, 15, rot = rot, gp = gpar(col=col, cex = GenomeGraphs:::getCex(gdObject)),
                      default.units = "native", just = c("center", "center"))
        }
    }
    popViewport(1)
}

setMethod("drawGD", signature("Gene"), .drawGene)
setMethod("drawGD", signature("GeneRegion"), .drawGene)

# =============================================================================
#			set single gene name to HGNC
# =============================================================================

setMethod("initialize", "Gene", function(.Object, ...){
    .Object <- callNextMethod()
    #.Object@ens <- getBM(c("structure_gene_stable_id","structure_transcript_stable_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","structure_transcript_chrom_strand", "structure_biotype"), filters = .Object@type, values=.Object@id, mart=.Object@biomart)
 .Object@ens <- getBM(c("external_gene_name","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand"), filters = .Object@type, values=.Object@id, mart=.Object@biomart)
    if(!is.null(.Object@ens)){
     .Object@ens <- cbind(.Object@ens, biotype=rep("protein_coding", length(.Object@ens[,1])))
   }
    if (is.null(.Object@ens)) {
        setPar(.Object, "size", 0)
    }
    .Object
})

makeGene <- function(id, type, biomart, dp = NULL){
 if(missing(id)) stop("Need to specify a gene identifier for creating a Gene")
  pt <- getClass("Gene")@prototype
 if (is.null(dp))
   dp <- pt@dp
 if(missing(type))
   type=pt@type
 new("Gene", id = id, type = type, biomart = biomart, dp = dp)
}


# ====================================================================================
#				FUNCTION
# ====================================================================================
createTranslocPlot = function(bamFile,CNfile,chrom,start,end,breakpoint,covFile=NULL,outFile=NULL,genome="grch37",plotChms="top",geneName=NULL,geneRot=90)
	{
	library(GenomeGraphs)
	# create region bed
	write.table(rbind(c(chrom,start,end,0,"+"),
			c(chrom,start,end,0,"-")),
		"region.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	# get coverage for region
	#if(is.null(covFile))
	#	{
	#	covFile=paste0("cov",chrom,".coverage")
	#	system(paste0("samtools depth ",bamFile," -b region.bed > ",covFile))
	#	}
	# region to plot
	region = read.table(
"region.bed",head=FALSE)
	# copy number
	print("copy number")
	if(rev(strsplit(CNfile,split="[.]")[[1]])[1]=="csv")
		{
		segs = read.csv(CNfile,head=FALSE)
		} else {
		segs = read.table(CNfile,head=FALSE,sep="\t")
		}
	if(ncol(segs)==6)
		{
		chromCol = 1
		segStartCol = 2
		segEndCol = 3
		segCNcol = 5
		segCNbCol = 6
		} else {
		chromCol = 2
		segStartCol = 3
		segEndCol = 4
		segCNcol = 7
		segCNbCol = 8
		}
	chromIndex = which(segs[,chromCol]==region[,1])
	seg = makeSegmentation(segs[chromIndex,segStartCol],
				segs[chromIndex,segEndCol],
				segs[chromIndex,segCNcol],
				dp = DisplayPars(color = "deepskyblue1", lwd=2,lty = 1))
	segA = makeSegmentation(segs[chromIndex,segStartCol],
				segs[chromIndex,segEndCol],
				segs[chromIndex,segCNcol]-segs[chromIndex,segCNbCol],
				dp = DisplayPars(color = "red", lwd=2,lty = 2))
	segB = makeSegmentation(segs[chromIndex,segStartCol],
				segs[chromIndex,segEndCol],
				segs[chromIndex,segCNbCol],
				dp = DisplayPars(color = "green", lwd=2,lty = 2))
	# coverage
	#cov = read.table(covFile,head=FALSE)
	#cop = makeGenericArray(intensity=as.matrix(cov[,3,drop=FALSE]),
	#		probeStart=cov[,2])
	# BAF
	print("BAF")
	if(rev(strsplit(bafFile,split="[.]")[[1]])[1]=="gz")
		{
		system(paste0("gunzip -c ",bafFile," | awk '$2==",region[,1],"&&$3>=",region[,2],"&&$3<=",region[,3],"' > BAF.txt"))
		} else {
		system(paste0("awk '$2==",region[,1],"&&$3>=",region[,2],"&&$3<=",region[,3],"' ",bafFile," > BAF.txt"))
		}
	baf = read.table("BAF.txt",head=FALSE)
	BAF = makeGenericArray(intensity=as.matrix(baf[,6,drop=FALSE]),probeStart=baf[,3], 
	dp = DisplayPars(size=3, color = "seagreen", type="dot"))
	# logR	
	print("logR")
	logR = makeGenericArray(intensity=as.matrix(baf[,4,drop=FALSE]),probeStart=baf[,3], 
	dp = DisplayPars(size=3, color = "orange", type="dot"))
	# raw copy number
	#print("copynumber")
	#cn = makeGenericArray(intensity=as.matrix(baf[,10,drop=FALSE]),
	#		probeStart=baf[,3],
	#		trackOverlay=list(seg,segA,segB), 
	#		dp = DisplayPars(size=5, color = "deepskyblue4", type="dot"))
	cn = makeGenericArray(intensity=as.matrix(c(0,max(segs[chromIndex,segCNcol])+1,rep(2,times=nrow(baf)-2))),
			probeStart=baf[,3],
			trackOverlay=list(seg,segA,segB), 
			dp = DisplayPars(size=5, color = "white", type="dot"))

	# axis etc
	print("ideogram")
	ideog = makeIdeogram(chromosome = region[,1])
	print("axis")
	genomeAxis <- makeGenomeAxis(add53 = TRUE, add35=TRUE)
	# biomart
	print("biomart")
	if(genome=="grch37")
	{
	  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
	                 host="grch37.ensembl.org", 
	                 path="/biomart/martservice",
	                 dataset="hsapiens_gene_ensembl")
	} else {
	mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
	}
	# begin plot params
	toPlot = list(ideog,
			genomeAxis)
	# genes
	print("genes")
	if(is.null(geneName))
		{
		genesplus <- makeGeneRegion(start=start, end=end,
			strand="+", chromosome=chrom, biomart=mart,
			dp = DisplayPars(plotId=TRUE, idRotation = geneRot,
			idColor = "red",size=2,cex=0.5))
		genesneg <- makeGeneRegion(start=start, end=end,
			strand="-", chromosome=chrom, biomart=mart,
			dp = DisplayPars(plotId=TRUE, idRotation = geneRot,
			idColor = "blue",size=2,cex=0.5))
		toPlot = append(toPlot,genesplus,genesneg)
		} else {
		gene = makeGene(geneName,type="hgnc_symbol",
			biomart=mart,dp=DisplayPars(plotId=TRUE,idRotation=geneRot,
						idColor="red",size=1,cex=0.5))
		toPlot = append(toPlot,gene)
		}
	# continue plot params
	toPlot = append(toPlot,
		    list(logR=logR,
		    BAF=BAF))
	names(toPlot)[[1]] = paste0("chms",chrom)
	if(plotChms!="top") toPlot = rev(toPlot)
	if(plotChms=="top")
		{
		region=c(2,6)
		if(!is.null(geneName)) region[2] = region[2]-1
		} else {
		region=c(1,5)
		if(!is.null(geneName)) region[2] = region[2]-1
		}
	# make plot
	print("plot")
	if(!is.null(outFile)) pdf(outFile)
		gdPlot(toPlot,
			overlay=c(makeRectangleOverlay(start=breakpoint,
						end=breakpoint,
						region=region,
						dp=DisplayPars(color="gray"))#,
				#makeRectangleOverlay(start=breakpoint,
				#		end=breakpoint,
				#		region=region2,
				#		dp=DisplayPars(color="gray"),
				#		coords="genomic")
))
	if(!is.null(outFile)) dev.off()
	}
