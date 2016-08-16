createSegFromConumee = function(dir,ignore=NULL)
	{
	files = list.files(dir)
	files = files[grep("-res",files)]
	if(!is.null(ignore))
		{
		index = unique(unlist(sapply(ignore,FUN=function(x) grep(x,files))))
		if(length(index)>0) files = files[-index]
		}
	out = NULL
	for(i in 1:length(files))
		{
		load(files[i])
		tmp = x@seg$summary
		tmp[,"seg.median"] = tmp[,"seg.median"]-x@bin$shift
		tmp[,1] = strsplit(files[i],split="-")[[1]][1]
		tmp[,2] = gsub("chr","",tmp[,2])
		tmp = tmp[,c(1:5,8)]
		out = rbind(out,tmp)
		}
	colnames(out) = c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
	write.table(out,file=paste0(dir,"/seg.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	}

segToMarkers=function(file,dir,manifestMarkersFile = "/media/grp11/CI_Pathology_Steele/methylation/manifests/EPICmarkers.txt")
	{
	columnNames = c("Marker Name", "Chromosome", "Marker Position")
	# read segment file 
	data = read.table(paste0(dir,'/',file),head=TRUE,sep='\t',as.is=TRUE)
	tmp1 = data[,c(2,3)]
	tmp2 = data[,c(2,4)]
	colnames(tmp2) = colnames(tmp1)
	out = rbind(tmp1,tmp2)
	out = out[order(out[,1],out[,2]),]
	out = out[!duplicated(out),]
	out = cbind(paste0('marker',1:nrow(out)),out)
	colnames(out) = columnNames
	# read manifest files
	manifest = read.table(manifestMarkersFile,sep="\t",head=TRUE,as.is=TRUE)
	colnames(manifest) = columnNames
	#dup = apply(out,MARGIN=1,FUN=function(x) any(manifest[,2]==x[2]&manifest[,3]==x[3]))
	#if(sum(dup)>0)
	#	{
	#	out = out[-which(dup),]
	#	}
	out = rbind(out,manifest)
	write.table(out,file=paste0(dir,'/markers.txt'),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
	}
