runMyChamp = function(dataDir,outDir,sampleSheet="Pillay EPIC samplesheet 280516.csv",studyInfo="studyInfo_2.txt",mychampFile="/home/chris/Dropbox/PostDoc/methylation/Rscripts/myCHAMP.R",min.width=2,nperm=10000,arrayType="EPIC",infoFactor=c(),controlGroup='C')
	{
	library(ChAMP)
	# load data
	data = champ.load(directory=dataDir,methValue="B",resultsDir=outDir,QCimages=TRUE,arraytype=arrayType)
	save(list="data",file=paste0(outDir,"/data.Rdata"))
	# normalise type 1 and 2 probes
	norm = champ.norm(beta=data$beta,rgSet=data$rgSet,pd=data$pd,mset=data$mset,sampleSheet=paste0(dataDir,"/",sampleSheet),resultsDir=outDir,methValue="B",fromIDAT=TRUE,norm="BMIQ",filter=TRUE,filterXY=TRUE,QCimages=TRUE,plotBMIQ=TRUE,arraytype=arrayType)
	save(list="norm",file=paste0(outDir,"/norm.Rdata"))
	# SVD
	champ.SVD(beta=norm$beta,rgSet=data$rgSet,detP=data$detP,pd=data$pd,sampleSheet=paste0(dataDir,"/",sampleSheet),studyInfo=TRUE,studyInfoFile=paste0(dataDir,"/",studyInfo),infoFactor=infoFactor,resultsDir=outDir)
	rm(norm)
	# CNA
	library(ChAMP)
	library(DNAcopy)
	source(mychampFile)
	mychamp.CNA(intensity=data$intensity,pd=data$pd,resultsDir=outDir,controlGroup=controlGroup,seg.min.width=min.width,seg.nperm=nperm,arraytype=arrayType)
	}




