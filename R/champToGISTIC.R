champToGistic=function(dir)
	{
	files = list.files(dir)
	files = files[grep(".txt",files)]
	out = NULL
	for(i in 1:length(files))
		{
		out = rbind(out,read.table(paste0(dir,"/",files[i]),head=TRUE))
		}
	out = out[,-ncol(out)]
	colnames(out) = c('Sample','Chromosome','Start Position','End Position','Num markers','Seg.CN')
	write.table(out,file=paste0(dir,'/seg.txt'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
	}

segToMarkers=function(file,dir)
	{
	data = read.table(paste0(dir,'/',file),head=TRUE,sep='\t')
	tmp1 = data[,c(2,3)]
	tmp2 = data[,c(2,4)]
	colnames(tmp2) = colnames(tmp1)
	out = rbind(tmp1,tmp2)
	out = out[order(out[,1],out[,2]),]
	out = out[!duplicated(out),]
	out = cbind(paste0('marker',1:nrow(out)),out)
	colnames(out) = c("Marker Name", "Chromosome", "Marker Position")
	write.table(out,file=paste0(dir,'/markers.txt'),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
	}



