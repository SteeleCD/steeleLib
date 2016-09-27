subsetVCF = function(vcfDat,chrom,pos)
	{
	index = mapply(FUN=function(c,p) which(vcfDat$CHROM==c&vcfDat$POS==p),c=chrom,p=pos)
	vcfDat = vcfDat[unlist(index),]
	return(vcfDat)
	}

getVcfData = function(normVcf,tumVcf)
	{
	tum = read.table(tumvCF,sep="\t",as.is=TRUE)
	norm = read.table(normVcf,sep="\t",as.is=TRUE)
	return(list(norm=norm,tum=tum))
	}

checkIdentity = function(dataT,dataN,SNPdata)
	{
	dataT = subsetVCF(dataT,SNPdata$Chromosome,SNPdata$Position..Hg19.)
	dataN = subsetVCF(dataN,SNPdata$Chromosome,SNPdata$Position..Hg19.)

	infoN = apply(dataN,MARGIN=1,FUN=function(x) paste0(x[1:3],collapse="-"))
	infoT = apply(dataT,MARGIN=1,FUN=function(x) paste0(x[1:3],collapse="-"))

	indexN = which(infoN%in%infoT)
	indexT = which(infoT%in%infoN)

	dataN = dataN[indexN,]
	dataT = dataT[indexT,]

	AFT = sapply(dataT$INFO,FUN=function(x) gsub("AF=","",strsplit(x,split="[;]")[[1]][2]))
	AFN = sapply(dataN$INFO,FUN=function(x) gsub("AF=","",strsplit(x,split="[;]")[[1]][2]))

	total = 2*length(AFT)
	matching = 2*length(which(AFT==AFN))
	index = which(AFT!=AFN)
	diffs = abs((2*AFT[index])-(2*AFN[index]))
	matching = matching+sum(diffs[which(diffs==1)])
	nonmatching = sum(diffs)
	matching/(matching+nonmatching)
	}

load.candidate.SNPs = function(file=NULL)
	{
	if(is.null(file))
		{
		tmpEnv = new.env()
		data(list="candidateSnps",package="steeleLib",envir=tmpEnv)
		return(tmpEnv[["candidateSnps"]])
		} else {
		data = read.table(file)
		return(data)
		}
	}


