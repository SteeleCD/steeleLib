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





