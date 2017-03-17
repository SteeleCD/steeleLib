# ===========================================================
# QOL function for rbinding many files
# ===========================================================

combineFiles = function(files,outFile,sep="\t",head=TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE)
	{
	out = NULL
	for(i in 1:length(files))
		{
		out = rbind(out,read.table(files[i],sep=sep,head=head))
		}
	write.table(out,file=outFile,col.names=col.names,row.names=row.names,sep=sep)
	}
