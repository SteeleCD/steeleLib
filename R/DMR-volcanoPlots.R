
plotVolcano = function(data,manifestProbes=NULL,file=NULL,PAR=list(mfrow=c(1,1),mar=c(5,4,4,2)),title="test",folder="/home/chris/Dropbox/PostDoc/methylation/plots/undiffSarc/DMP/volcano/",MAIN="test")
	{
	base = paste0(folder,title,"/")
	if(!is.null(file)) pdf(paste0(base,title,file))
	if(!is.null(PAR)) do.call(par,PAR)
	index = which(data$adj.P.Val<0.05)
	Xvals = data$logFC[index]
	Yvals = -log(data$adj.P.Val[index])
	probes = paste0(data$probeID[index])
	XLIM = range(Xvals)
	YLIM = range(Yvals)
	if(!is.null(manifestProbes))
		{
		manifestProbes = paste0(manifestProbes)
		Xvals = Xvals[which(probes%in%manifestProbes)]
		Yvals = Yvals[which(probes%in%manifestProbes)]
		probes = probes[which(probes%in%manifestProbes)]
		}
	plot(Xvals,Yvals,col=ifelse(abs(Xvals)<0.3,"black","red"),ylab="-log(P)",xlab="log(fold change)",xlim=XLIM,ylim=YLIM,main=MAIN)
	if(!is.null(file)) dev.off()
	return(list(hyper=probes[which(Xvals>0.3)],hypo=probes[which(Xvals<c(-0.3))]))
	}
