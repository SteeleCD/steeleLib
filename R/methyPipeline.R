methPipeline = function(dataDir,outDir=getwd(),nCells=5)
{
# load in sample sheet
sheet = read.metharray.sheet(dataDir)
# read idats based on sample sheet
RGset = read.metharray.exp(targets=sheet)
# detection p values
qualCut = 0.01
detP = detectionP(RGset)
lowQual = detP>qualCut
# fraction of failed positions per sample
colMeans(lowQual) 
probesToKeep = rownames(detP)[which(rowSums(detP>=qualCut)<=0)]
# functional normalisation
funnorm = preprocessFunnorm(RGset)
# drop SNP probes
funnorm = dropLociWithSnps(funnorm)
save(list="funnorm",file=paste0(outDir,"funnorm-droppedSnps.Rdata"))
# get betas
betas = getBeta(funnorm)
manifest = getManifest(RGset)
annotation = getAnnotation(funnorm)
# remove sex probes
sexProbes = rownames(annotation)[which(annotation[,"chr"]%in%c("chrX","chrY"))]
betas = betas[-which(rownames(betas)%in%sexProbes),]
# keep high quality probes
betas = betas[which(rownames(betas)%in%probesToKeep),]
# save betas
colnames(betas) = sheet$Sample_Name
save(list="betas",file=paste0(outDir,"betas-funnorm-dropSex-dropLowQ.Rdata"))

# create refactor file
tmp = rbind(rownames(betas),betas)
colnames(tmp)[1] = "ID"
write.table(tmp,file=paste0(outDir,"refactor_input.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)

# run refactor
refactorOutput <- refactor(paste0(outDir,"refactor_input.txt"), nCells, out = "demo_refactor")

# make factors
batch = as.factor(sheet$Sentrix_ID)
covModel = model.matrix(~refactorOutput$refactor_components)
# run combat
combatMs = ComBat(dat=logit(undiffBetas),batch=batch,mod=covModel)
# save combat Ms
save(list="combatMs",file=paste0(outDir,"combatMs.Rdata"))
}