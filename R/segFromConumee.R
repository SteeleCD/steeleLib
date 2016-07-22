createSegFromConumee = function(dir)
	{
	files = list.files(dir)
	files = files[grep("-res",files)]
	out = NULL
	for(i in 1:length(files))
		{
		load(files[i])
		tmp = x@seg$summary
		tmp[,1] = strsplit(files[i],split="-")[[1]][1]
		tmp[,2] = gsub("chr","",tmp[,2])
		out = rbind(out,tmp)
		}
	write.table(out,file=paste0(dir,"/seg.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	}
