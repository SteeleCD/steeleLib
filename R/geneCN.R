# function to determine whether a gene is amplified or deleted based on log2(CN/ploidy)
getGeneAmpDel = function(seg,			# seg file
			chr,			# gene chromosome
			start,			# gene start
			end,			# gene end
			type="amplified",	# "amplified" or "deleted"
			threshold=1.3,		# CN threshold
			CNcol=15,		# log2(CN/ploidy) seg column
			chromCol=3,		# chromosome seg column
			startCol=4,		# start seg column
			endCol=5		# end seg column
			)
	{
	seg.info = paste0("chr",seg[,chromCol],":",seg[,startCol],"-",seg[,endCol])
	seg.gr = as(seg.info,"GRanges")
	gene.gr = as(paste0(chr,":",start,"-",end),"GRanges")
	overlap = findOverlaps(seg.gr,gene.gr)
	seg.overlap = seg[overlap@from,]
	if(type=="amplified")
		{
		return(min(seg.overlap[,CNcol])>threshold)
		} else {
		return(min(seg.overlap[,CNcol])<c(-threshold))
		}
	}

# function to determine whether a gene has LOH
getGeneLOH = function(seg,			# seg file
			chr,			# gene chromosome
			start,			# gene start
			end,			# gene end
			totCNcol=12,		# total CN column
			aCNcol=13,		# A allele CN column
			bCNcol=14,		# B allele CN column
			chromCol=3,		# chromosome seg column
			startCol=4,		# start seg column
			endCol=5		# end seg column
			)
	{
	CNcols=c(totCNcol,aCNcol,bCNcol)
	seg.info = paste0("chr",seg[,chromCol],":",seg[,startCol],"-",seg[,endCol])
	seg.gr = as(seg.info,"GRanges")
	gene.gr = as(paste0(chr,":",start,"-",end),"GRanges")
	overlap = findOverlaps(seg.gr,gene.gr)
	seg.overlap = seg[overlap@from,]
	bools = apply(seg.overlap[,CNcols],MARGIN=1,FUN=function(x) any(x[2:3]==0)&x[1]!=0)
	return(any(bools))
	}

# function to determine whether a gene has homozygous deletion
getGeneHomDel = function(seg,			# seg file
			chr,			# gene chromosome
			start,			# gene start
			end,			# gene end
			CNcol=15,		# log2(CN/ploidy) seg column
			chromCol=3,		# chromosome seg column
			startCol=4,		# start seg column
			endCol=5		# end seg column
			)
	{
	seg.info = paste0("chr",seg[,chromCol],":",seg[,startCol],"-",seg[,endCol])
	seg.gr = as(seg.info,"GRanges")
	gene.gr = as(paste0(chr,":",start,"-",end),"GRanges")
	overlap = findOverlaps(seg.gr,gene.gr)
	seg.overlap = seg[overlap@from,]
	return(any(seg.overlap[,CNcol]==0))
	}

# function to obtain CN at a gene
getGeneCN = function(seg,			# seg file
			chr,			# gene chromosome
			start,			# gene start
			end,			# gene end
			totCNcol=12,		# total CN column
			aCNcol=13,		# A allele CN column
			bCNcol=14,		# B allele CN column
			chromCol=3,		# chromosome seg column
			startCol=4,		# start seg column
			endCol=5		# end seg column
			)
	{
	CNcols=c(totCNcol,aCNcol,bCNcol)
	seg.info = paste0("chr",seg[,chromCol],":",seg[,startCol],"-",seg[,endCol])
	seg.gr = as(seg.info,"GRanges")
	gene.gr = as(paste0(chr,":",start,"-",end),"GRanges")
	overlap = findOverlaps(seg.gr,gene.gr)
	seg.overlap = seg[overlap@from,]
	return(seg.overlap)
	}

# function to get gene information
getGeneInfo = function(genesIn)
	{
	# get gene names
	library(org.Hs.eg.db)
	unmapped = org.Hs.eg.db::org.Hs.egSYMBOL
	mapped = mappedkeys(unmapped)
	genes = unlist(as.list(unmapped[mapped]))
	# get transcripts
	library(Homo.sapiens)
	hs = Homo.sapiens
	tr = transcriptsBy(hs, by="gene")
	# change gene names to hugo gene names
	hugo = sapply(names(tr),FUN=function(x) genes[which(names(genes)==x)])
	names(tr) = hugo
	# change from transcript to gene
	gene.gr <- reduce(tr) # ISA GenomicRange
	gene.info<-as(gene.gr,'data.frame') 
	# subset to genes of interest
	theseGenes = which(gene.info$group_name%in%genesIn)
	gene.info = gene.info[theseGenes,]
	return(gene.info)
	}
