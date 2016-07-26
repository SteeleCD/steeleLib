methPipeline = function(dataDir,outDir=getwd(),nCells=5,legendloc="topleft",plotcolours=c("green","red","black","blue","cyan"))
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
tmp = cbind(rownames(betas),betas)
colnames(tmp)[1] = "ID"
write.table(tmp,file=paste0(outDir,"refactor_input.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)
rm(tmp)
# run refactor to estimate cell proportions in each sample
refactorOutput <- refactor(paste0(outDir,"refactor_input.txt"), nCells, out = paste0(outDir,"refactor_out.Rdata"))
# make factors for combat
batch = as.factor(sheet$Slide)
covModel = model.matrix(~refactorOutput$refactor_components)
# run combat to normalise based on batch and cell proportions
combatMs = ComBat(dat=logit(betas),batch=batch,mod=covModel)
# save combat Ms
save(list="combatMs",file=paste0(outDir,"combatMs.Rdata"))
# pca
pca = runPCA(data=combatMs,fileOut=paste0(outDir,"pca.Rdata"))
plotPCA(pca,fileName="pca.pdf",outDir=outDir,groups=sheet$Pool_ID,arrows=FALSE,legendloc=legendloc,colours=plotcolours)
# hclust
groups = sapply(levels(as.factor(sheet$Pool_ID)),FUN=function(x) sheet$Sample_Name[which(sheet$Pool_ID==x)])
hClust = runHclust(pca$x,groups=groups,file=paste0(outDir,"hclust.pdf"),plotcolours=colours)
save(list="hClust",file=paste0(outDir,"hclust.Rdata"))
}
